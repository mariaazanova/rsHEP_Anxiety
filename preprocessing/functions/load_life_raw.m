function [EEG, ECG] = load_life_raw(subjid, edf_dir, ecg_label, eog_labels, chanlocs)
f = fullfile(edf_dir, [subjid, '.edf']);
data = edfread(f);
info = edfinfo(f);
EEG = eeg_emptyset;
EEG.data = single(cell2mat(cellfun(@(c) vertcat(data.(c){:}), data.Properties.VariableNames, 'UniformOutput', false))');
EEG.srate = info.NumSamples(1) / seconds(info.DataRecordDuration);
EEG.chanlocs = struct('labels', cellstr(info.SignalLabels)');
EEG = eeg_checkset(EEG);

% ECG: keep only the ECG channel
ECG = pop_select(EEG, 'channel', find(strcmp({EEG.chanlocs.labels}, ecg_label)));
ECG.chanlocs(1).labels = 'ECG';

% EEG: drop ECG and EOG channels, look up channel locations
EEG = pop_select(EEG, 'nochannel', find(ismember({EEG.chanlocs.labels}, [{ecg_label}, eog_labels])));
% one recording labels FT9/FT10 as T1/T2
labels = {EEG.chanlocs.labels};
if any(strcmp(labels, 'T1')), EEG.chanlocs(strcmp(labels, 'T1')).labels = 'FT9';  end
if any(strcmp(labels, 'T2')), EEG.chanlocs(strcmp(labels, 'T2')).labels = 'FT10'; end
EEG = pop_chanedit(EEG, 'lookup', chanlocs);
end
