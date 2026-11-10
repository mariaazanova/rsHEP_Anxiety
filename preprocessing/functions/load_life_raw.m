function [EEG, ECG] = load_life_raw(subjid, raw_dirs, ecg_label, eog_labels, chanlocs)
f = [];
for d = 1:length(raw_dirs)
    f = [f; dir(fullfile(raw_dirs{d}, '**', [subjid, '*.vhdr']))];
end

% prefer the _2 recording if they exist (correspond to restarts of the recording)
% then <ID>_Ruhe.vhdr, then whatever exists (typos)
second = ~cellfun(@isempty, regexpi({f.name}, ['^', subjid, '(_Ruhe_2|_2_Ruhe)\.vhdr$'], 'once'));
regular = strcmpi({f.name}, [subjid, '_Ruhe.vhdr']);
if any(second)
    f = f(find(second, 1));
elseif any(regular)
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
drop = ismember({EEG.chanlocs.labels}, [{ecg_label}, eog_labels, {'nase'}]);
EEG = pop_select(EEG, 'nochannel', find(drop));

% some ESF recordings label FT9/FT10 as T1/T2
labels = {EEG.chanlocs.labels};
if any(strcmp(labels, 'T1')), EEG.chanlocs(strcmp(labels, 'T1')).labels = 'FT9';  end
if any(strcmp(labels, 'T2')), EEG.chanlocs(strcmp(labels, 'T2')).labels = 'FT10'; end
EEG = pop_chanedit(EEG, 'lookup', chanlocs);
end
