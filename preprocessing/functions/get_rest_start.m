function [t0, t1, src] = get_rest_start(EEG, subjid, christian_dirs, rest_min)
% GET_REST_START  Start and end (s, in EEG's time base) of the resting-state block.
%   The block runs from marker S99 to S100.
%   src: 'S99'       - S99 in the recording; without S100, t1 = t0 + rest_min
%        'S100'      - only S100 in the recording: t0 = S100 - rest_min
%        'christian' - no marker in the recording: S99 position from Christian's marker files
%                      (<ID>_Ruhe_Startmarker_S99 in one of christian_dirs; header-only files are ignored)
%        'none'      - no start marker anywhere: t0 = t1 = NaN, the caller keeps the whole recording
%   t0 can be negative, or t1 beyond the recording; the caller keeps the part that was recorded.
types = cellfun(@(t) strrep(char(t), ' ', ''), {EEG.event.type}, 'UniformOutput', false);
lat = [EEG.event.latency];
i99 = find(strcmp(types, 'S99'));
i100 = find(strcmp(types, 'S100'));
if ~isempty(i99)
    t0 = (lat(i99) - 1) / EEG.srate;
    t1 = t0 + rest_min * 60;
    if ~isempty(i100), t1 = (lat(i100) - 1) / EEG.srate; end
    src = 'S99';
elseif ~isempty(i100)
    t1 = (lat(i100) - 1) / EEG.srate;
    t0 = t1 - rest_min * 60;
    src = 'S100';
else
    t0 = christian_start(subjid, christian_dirs);
    t1 = t0 + rest_min * 60;
    src = 'christian';
    if isnan(t0), src = 'none'; end
end
end

function t0 = christian_start(subjid, dirs)
t0 = NaN;
for d = 1:numel(dirs)
    f = dir(fullfile(dirs{d}, [subjid, '*']));
    f = f(~[f.isdir] & ~cellfun(@isempty, regexpi({f.name}, ['^', subjid, '_+Ruhe_Startmarker_S99$'], 'once')));
    for k = 1:numel(f)
        txt = fileread(fullfile(f(k).folder, f(k).name));
        p = regexp(txt, 'S 99,\s*(\d+)', 'tokens', 'once');
        if ~isempty(p)          % header-only files hold no S99 marker
            r = regexp(txt, 'Sampling rate:\s*(\d+)', 'tokens', 'once');
            t0 = str2double(p{1}) / str2double(r{1});
            return
        end
    end
end
end
