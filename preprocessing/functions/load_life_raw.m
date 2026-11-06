function [EEG, ECG] = load_life_raw(subjid, raw_dirs, ecg_label, eog_labels, chanlocs)
f = [];
for d = 1:length(raw_dirs)
    f = [f; dir(fullfile(raw_dirs{d}, '**', [subjid, '*.vhdr']))];
end

%some subjects have typos or only _2 recordings, ideally we take the first one,
% or whatever else is available. 
regular = strcmpi({f.name}, [subjid, '_Ruhe.vhdr']);
if any(regular)
    f = f(find(regular, 1));
else
    f = f(1);
end
EEG = pop_loadbv(f.folder, f.name);

% ECG: keep only the ECG channel
ecg_idx = find(strcmp({EEG.chanlocs.labels}, ecg_label));
ECG = pop_select(EEG, 'channel', ecg_idx);
ECG.chanlocs(1).labels = 'ECG';

% EEG: drop ECG and EOG channels, look up channel locations
drop = ismember({EEG.chanlocs.labels}, [{ecg_label}, eog_labels]);
EEG = pop_select(EEG, 'nochannel', find(drop));
EEG = pop_chanedit(EEG, 'lookup', chanlocs);
end
