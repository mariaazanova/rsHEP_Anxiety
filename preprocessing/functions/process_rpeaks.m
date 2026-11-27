function [qrs, is_ext, info] = process_rpeaks(ecg, qrs, fs, epoch_s, n_ext_cers, hr_cers, p)
% PROCESS_RPEAKS  From HEPLAB R-peaks to the R-peaks used for the HEP, and which of them are extrasystolic.
%   ecg        : ECG signal (vector, any polarity);  qrs : HEPLAB R-peak sample indices;  fs : sampling rate
%   epoch_s    : HEP epoch [start end] in s;  n_ext_cers, hr_cers : CER-S V+S count and heart rate (NaN if none)
%   p          : rpeak_params from config.m
%   1) Align: HEPLAB marks the largest deflection near each QRS, which in biphasic
%      complexes can jump between R and S wave. Each marker is moved to the same point of the participant's
%      median QRS template (best correlation within +-align_shift_ms); beats matching it less than
%      align_min_corr (extrasystoles, noise) keep their marker. Markers within 100 ms are one QRS: the
%      better match is kept.
%   2) Unreliable ECG: 10-s windows where HEPLAB and an independent R-peak
%      detector agree on fewer than min_bsqi of the beats (bSQI, Li, Mark & Clifford 2008), or the ECG is lost
%      (< min_beats R-peaks), flat or clipped. R-peaks whose epoch overlaps such a stretch are not used.
%   3) Extrasystoles (We have CER-S counts only). CER-S is trusted if its heart rate is within
%      max_hr_dev_cers of HEPLAB's on clean ECG stretches; then 0 CER-S extrasystoles means none. Otherwise,
%      or if CER-S found any, a beat is extrasystolic if premature (RR < ratio_thr x running reference),
%      abnormally shaped (corr < corr_thr with the normal-beat template, +-ext_shift_ms), or next to a gap
%      (RR > gap_thr x reference: a missed beat). The beat after each extrasystole is flagged too (according to RegRep).
%   qrs : R-peaks used;  is_ext : logical per qrs;  info : intermediate results for extrasystoles.mat

info = struct('qrs_heplab', qrs(:), 'n_ext_cers', n_ext_cers, 'hr_cers', hr_cers);
x = double(ecg(:));
[b, a] = butter(2, [0.5 40] / (fs/2), 'bandpass');
xf = filtfilt(b, a, x);                           % 0.5-40 Hz: QRS windows and detector peak positions

% 1) align
[qrs, info.align] = align(xf, qrs(:), fs, p);

% 2) unreliable ECG
[segs, info.quality] = bad_stretches(x, xf, qrs, fs, p);
ep = round(epoch_s * fs);
in_bad = false(numel(qrs), 1);
for k = 1:size(segs, 1)
    in_bad = in_bad | (qrs + ep(2) >= segs(k, 1) & qrs + ep(1) <= segs(k, 2));
end
qrs = qrs(~in_bad);
info.bad_segs = segs; info.n_dropped_bad_ecg = sum(in_bad);

% 3) extrasystoles
info.hr_heplab = numel(qrs) / ((numel(x) - sum(diff(segs, 1, 2) + 1)) / fs / 60);
info.cers_ok = abs(info.hr_heplab / hr_cers - 1) <= p.max_hr_dev_cers;   % false without CER-S result (NaN)
if info.cers_ok && n_ext_cers == 0
    is_ext = false(numel(qrs), 1);
    info.rule = [];
else
    [is_ext, info.rule] = extrasystoles(xf, qrs, fs, p);
    is_ext = is_ext | [false; is_ext(1:end-1)];
end
info.qrs = qrs; info.is_ext = is_ext;
end

%% 1) alignment
function [qrs_out, info] = align(xf, qrs, fs, p)
L = round(p.align_shift_ms / 1000 * fs);
[W, inside] = qrs_windows(xf, qrs, fs, L);
pre = round(0.080 * fs); n = size(W, 2) - 2 * L;
tmpl = median(W(inside, L + 1 : L + n), 1);       % median of all beats, then twice of the well-matching ones
for it = 1:2
    [c, ~, A] = best_corr(W, tmpl, L);
    good = c > p.align_min_corr;
    tmpl = median(A(good, :), 1);
end
[c, lag] = best_corr(W, tmpl, L);
[~, fid] = max(abs(tmpl));                         % marker = template's largest absolute deflection
move = c >= p.align_min_corr;
shift = zeros(numel(qrs), 1);
shift(move) = lag(move) + fid - pre - 1;
qrs_out = qrs + shift;
kept = true(numel(qrs), 1);                        % merge double detections (< 100 ms apart)
cc = c; cc(isnan(cc)) = -Inf;
while true
    idx = find(kept);
    d = find(diff(qrs_out(idx)) < round(0.100 * fs), 1);
    if isempty(d), break, end
    pair = idx(d:d + 1);
    [~, worse] = min(cc(pair));
    kept(pair(worse)) = false;
end
qrs_out = qrs_out(kept);
info = struct('shift', shift, 'corr', c, 'kept', kept, 'n_merged', sum(~kept));
end

%% 2) unreliable ECG
function [segs, info] = bad_stretches(x, xf, qrs, fs, p)
det = energy_detect(x, xf, fs);
tol = floor(p.bsqi_tol_s * fs);
w = round(p.win_s * fs);
starts = 1:round(w / 2):max(1, numel(x) - w + 1);
sd = arrayfun(@(s) std(x(s:s + w - 1)), starts(:));
xmin = min(x); xmax = max(x);
reason = false(numel(starts), 4);                  % lost, low bSQI, flat, clipped
bsqi = nan(numel(starts), 1);
for k = 1:numel(starts)
    s = starts(k); e = s + w - 1; seg = x(s:e);
    a = qrs(qrs >= s & qrs <= e); b = det(det >= s & det <= e);
    bsqi(k) = beat_f1(a, b, tol);
    reason(k, :) = [numel(a) < p.min_beats, bsqi(k) < p.min_bsqi, sd(k) < p.flat_rel * median(sd), ...
                    mean(seg <= xmin | seg >= xmax) > p.clip_frac];
end
mask = false(numel(x), 1);
for k = find(any(reason, 2))'
    mask(starts(k):starts(k) + w - 1) = true;
end
d = diff([false; mask; false]);
segs = [find(d == 1), find(d == -1) - 1];
info = struct('bsqi', bsqi, 'n_lost', sum(reason(:, 1)), 'n_low_bsqi', sum(reason(:, 2)), 'n_flat', sum(reason(:, 3)), ...
              'n_clipped', sum(reason(:, 4)), 'bad_fraction', mean(mask));
end

function f = beat_f1(a, b, tol)
% F1 agreement of two sets of beat positions
if isempty(a) && isempty(b), f = 1; return, end
if isempty(a) || isempty(b), f = 0; return, end
m = sum(min(abs(a(:) - b(:)'), [], 2) <= tol);
f = 2 * m / (numel(a) + numel(b));
end

function pk = energy_detect(x, xf, fs)
% Pan-Tompkins-like detector: 5-20 Hz band-pass, derivative, squaring, 150-ms integration; peaks above
% 0.3 x the local (5-s) 95th percentile, 250-ms refractory (higher peak wins), moved to the largest
% absolute deflection of the 0.5-40 Hz ECG within +-100 ms
[b, a] = butter(2, [5 20] / (fs/2), 'bandpass'); y = filtfilt(b, a, x);
n = floor(0.15 * fs); e = conv([0; diff(y)].^2, ones(n, 1) / n, 'same');
dec = 5; ed = e(1:dec:end); half = floor(fs * 5 / dec / 2);
lvl = zeros(size(ed));
for i = 1:numel(ed)
    lvl(i) = prctile(ed(max(1, i - half):min(numel(ed), i + half)), 95);
end
lvl = repelem(lvl, dec); lvl = lvl(1:numel(e));
e(e < 0.3 * lvl) = 0;
[~, pk] = findpeaks(e, 'MinPeakDistance', floor(0.25 * fs));
h = floor(0.1 * fs);
for i = 1:numel(pk)
    s = max(1, pk(i) - h); [~, j] = max(abs(xf(s:min(numel(xf), pk(i) + h - 1)))); pk(i) = s + j - 1;
end
end

%% 3) extrasystole rule 
function [is_ext, info] = extrasystoles(xf, qrs, fs, p)
ratio = [NaN; ref_ratio(diff(qrs) / fs * 1000, p.accept_thr, p.ref_win, p.max_rejected)];
next_ratio = [ratio(2:end); NaN];
L = round(p.ext_shift_ms / 1000 * fs);
[W, inside] = qrs_windows(xf, qrs, fs, L);
n = size(W, 2) - 2 * L;
% template: median of on-time beats not followed by a pause
normal = inside & ratio > 0.9 & ratio < 1.1 & next_ratio > 0.9 & next_ratio < 1.15;
tmpl = median(W(normal, L + 1 : L + n), 1);
c = best_corr(W, tmpl, L);
by_ibi = ratio < p.ratio_thr;                      % NaN compares false
by_shape = c < p.corr_thr;
gap = ratio > p.gap_thr;
by_gap = gap | [gap(2:end); false];
is_ext = by_ibi | by_shape | by_gap;
info = struct('ratio', ratio, 'corr', c, 'by_ibi', by_ibi, 'by_shape', by_shape, 'by_gap', by_gap, 'n_template_beats', sum(normal));
end

function ratio = ref_ratio(rr, accept_thr, win, max_rejected)
% RR / running median of recently accepted RR intervals; premature intervals and the pause after them
% are not accepted; after max_rejected rejections in a row the reference is rebuilt from a stable stretch
ratio = nan(size(rr));
acc = rr(1:min(win, numel(rr)));
rejected = 0; after_premature = false;
for i = 1:numel(rr)
    ref = median(acc(max(1, end - win + 1):end));
    ratio(i) = rr(i) / ref;
    if rr(i) < accept_thr * ref || rr(i) > 1.5 * ref
        rejected = rejected + 1;
        after_premature = rr(i) < accept_thr * ref;
    elseif after_premature
        after_premature = false;
        rejected = rejected + 1;
    else
        acc(end + 1) = rr(i); 
        rejected = 0;
    end
    if rejected >= max_rejected
        recent = rr(max(1, i - win + 1):i);
        if max(recent) - min(recent) < 0.2 * median(recent)
            acc = recent;
            rejected = 0;
        end
    end
end
end

%% shared: QRS windows and shifted template correlation
function [W, inside] = qrs_windows(xf, qrs, fs, L)
% QRS segments (-80..+120 ms) of the filtered ECG 
pre = round(0.080 * fs); post = round(0.120 * fs);
inside = qrs > pre + L & qrs <= numel(xf) - post - L;
W = nan(numel(qrs), pre + post + 2 * L);
for k = find(inside)'
    W(k, :) = xf(qrs(k) - pre - L : qrs(k) + post + L - 1);
end
end

function [c, lag, A] = best_corr(W, tmpl, L)
% highest correlation of each window with the template over shifts -L..+L, the shift, and the shifted windows
n = size(tmpl, 2);
t = tmpl - mean(tmpl);
C = nan(size(W, 1), 2 * L + 1);
for s = -L:L
    B = W(:, L + 1 + s : L + s + n);
    B = B - mean(B, 2);
    C(:, s + L + 1) = (B * t') ./ (sqrt(sum(B.^2, 2)) * norm(t));
end
[c, j] = max(C, [], 2);
c(all(isnan(C), 2)) = NaN;
lag = j - L - 1;
if nargout > 2
    A = nan(size(W, 1), n);
    for k = find(~isnan(c))'
        A(k, :) = W(k, L + 1 + lag(k) : L + lag(k) + n);
    end
end
end
