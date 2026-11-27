%% Alpha power per channel for the four ICA cleaning levels (protocol 2.3.5)
% Average power in 8-12 Hz from Welch PSD estimates (4 s Hanning windows, 50% overlap), per channel, computed
% four times: no ICA cleaning, no CFA ICs removed, 1-2 CFA ICs removed, 3-4 CFA ICs removed.
%
% a1-a4 only save epoched data (-200 to 600 ms), which cannot give 4 s windows. This script therefore rebuilds
% the continuous preprocessed recording exactly as a1 does (low-pass, notch, resample, high-pass, removal and
% interpolation of the channels a1 removed, average reference) and applies the ICA weights and the component
% selections that a2/a4 stored in the cleaned datasets.
%
% CHOICE not fixed by the protocol: alpha power is computed for one high-pass setting only (alpha_hp in config.m);
% the protocol states that alpha power depends on the ICA level only.
%
% Output (derivative_dir): alpha_ic_<level>_all_channels.csv - subject_id, one column per channel, roi_alpha
% (ROI = same electrodes as the HEP ROI). Units: uV^2/Hz.
% Test runs: set subj_limit before calling; output then goes to <deriv_root>/test/.

run('config.m');   % alpha_band, welch_window_s, welch_overlap, alpha_hp

files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); 
out_dir = derivative_dir;
if exist('subj_limit', 'var'), files = files(1:min(subj_limit, end)); out_dir = [deriv_root, 'test/']; end
create_dirs({out_dir});

n_levels = length(ic_cutoffs);
alpha_local = cell(length(files), 1);

parfor i = 1:length(files)
    try
        subjid = files(i).name;
        prep_eegpath = fullfile(files(i).folder, files(i).name, 'prep');
        clean_path = fullfile(files(i).folder, files(i).name, 'clean');

        % epoched a1 output: tells which channels a1 removed for this high-pass setting
        prep_id = sprintf('%s_hp_%0.1f.set', subjid, hp_settings(alpha_hp));
        PREP = pop_loadset('filename', prep_id, 'filepath', prep_eegpath, 'loadmode', 'info');

        % continuous data, same steps as a1
        EEG = load_life_raw(subjid, edf_dir, ecg_label, eog_labels, orig_chanlocs);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0, 'hicutoff', lp_settings, 'filtorder', lp_filter_order, 'plotfreqz', 0);
        EEG = pop_eegfiltnew(EEG, 'locutoff', notch_band(1), 'hicutoff', notch_band(2), 'revfilt', 1, 'filtorder', notch_filter_order, 'plotfreqz', 0);
        EEG = pop_resample(EEG, fs);
        chanlocs = EEG.chanlocs;
        EEG = pop_eegfiltnew(EEG, 'locutoff', hp_settings(alpha_hp), 'hicutoff', 0, 'filtorder', hp_filter_order(alpha_hp), 'plotfreqz', 0);
        if ~isempty(PREP.reject.removed_channels)
            EEG = pop_select(EEG, 'nochannel', PREP.reject.removed_channels);
        end
        EEG = pop_interp(EEG, chanlocs, 'spherical');
        EEG = pop_reref(EEG, []);

        % original ICA decomposition (a2). The cleaned datasets cannot be used for the weights: pop_subcomp in a4
        % drops the removed components from them, while etc.ic_remove still indexes the original decomposition.
        ica_id = sprintf('%s_hp_%0.1f_ica.set', subjid, hp_settings(alpha_hp));
        ICA = pop_loadset('filename', ica_id, 'filepath', fullfile(files(i).folder, files(i).name, 'ica'), 'loadmode', 'info');
        % ICA weight columns follow the a2 channel order
        if ~isequal({EEG.chanlocs.labels}, {ICA.chanlocs.labels}), error('Channel order differs from the ICA dataset'); end
        [has_chan, chan_order] = ismember(ref_labels, {EEG.chanlocs.labels});
        if ~all(has_chan), error('Missing channels: %s', strjoin(ref_labels(~has_chan), ', ')); end

        win = round(welch_window_s * EEG.srate);
        a = nan(n_levels, EEG.nbchan);
        for ic = 1:n_levels
            if isnan(ic_cutoffs(ic))
                X = EEG;                                       % no ICA cleaning
            else
                clean_id = sprintf('%s_hp_%0.1f_ic_%0.1f.set', subjid, hp_settings(alpha_hp), ic_cutoffs(ic));
                C = pop_loadset('filename', clean_id, 'filepath', clean_path, 'loadmode', 'info');
                r = C.etc.ic_remove;
                reject = unique([r.heart(:)', r.muscle(:)', r.eye(:)', r.LineNoise(:)', r.ChannelNoise(:)', r.other(:)']);
                X = EEG;
                X.icaweights = ICA.icaweights; X.icasphere = ICA.icasphere; X.icawinv = ICA.icawinv; X.icachansind = ICA.icachansind;
                if size(C.icaweights, 1) ~= size(ICA.icaweights, 1) - numel(reject)   % stale ica/ vs clean/ files
                    error('%s: a4 removed %d components, but ic_remove lists %d', clean_id, size(ICA.icaweights, 1) - size(C.icaweights, 1), numel(reject));
                end
                if ~isempty(reject), X = pop_subcomp(X, reject); end
            end
            [pxx, f] = pwelch(double(X.data)', hanning(win), round(welch_overlap * win), win, X.srate);   % columns = channels
            a(ic, :) = mean(pxx(f >= alpha_band(1) & f <= alpha_band(2), :), 1);
        end
        alpha_local{i} = a(:, chan_order);   % columns in ref_labels order
        if isfile([save_error, subjid, '_alpha_error_log.txt']), delete([save_error, subjid, '_alpha_error_log.txt']); end
    catch ME
        % failed subjects stay empty -> error log per subject
        fileID = fopen([save_error, files(i).name, '_alpha_error_log.txt'], 'w');
        fprintf(fileID, '%s\n', getReport(ME));
        fclose(fileID);
    end
end

valid = ~cellfun(@isempty, alpha_local);
[~, roi_idx] = ismember(electrodes, ref_labels);

subject_ids = {files(valid).name}';
for ic = 1:n_levels
    A = cell2mat(cellfun(@(a) a(ic, :), alpha_local(valid), 'UniformOutput', false));
    T = [table(subject_ids, 'VariableNames', {'subject_id'}), array2table(A, 'VariableNames', ref_labels), ...
         table(mean(A(:, roi_idx), 2), 'VariableNames', {'roi_alpha'})];
    writetable(T, sprintf('%salpha_ic_%0.1f_all_channels.csv', out_dir, ic_cutoffs(ic)));
end
fprintf('Alpha power: %d subjects written, %d failed\n', sum(valid), sum(~valid));
