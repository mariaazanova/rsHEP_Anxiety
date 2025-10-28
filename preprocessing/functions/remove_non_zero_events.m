function EEG = remove_non_zero_events_ECG(EEG)
%Cleans up the EEG structure by keeping only
% events at latency zero and removing duplicates and empty epochs.

for kk = 1:size(EEG.epoch,2)

    if isfield(EEG, 'epoch') && numel(EEG.epoch(kk).event) > 1

        keep_evnt_idx = cell2mat(EEG.epoch(kk).eventlatency) == 0;

        event_fields = {'event', 'eventbepoch', 'eventbini', 'eventbinlabel', 'eventbvmknum', 'eventbvtime', ...
            'eventchannel', 'eventcodelabel', 'eventduration', 'eventenable', 'eventflag', ...
            'eventitem', 'eventlatency', 'eventtype', 'eventvisible'};

        for i = 1:length(event_fields)
            field_name = event_fields{i};
            if isfield(EEG.epoch(kk), field_name)
                EEG.epoch(kk).(field_name) = EEG.epoch(kk).(field_name)(keep_evnt_idx);
            end
        end

    end

end

%events to keep
EEG.event = EEG.event([EEG.epoch.event]);
EEG = eeg_checkset(EEG,  'eventconsistency' );

% remove duplicate events (since events are multiplied during epoching, remove the arbitrary ones and keep one event per epoch)
if isfield(EEG.event, 'item')
    [~, ia, ~] = unique([EEG.event.item]);
    EEG.event = EEG.event(ia);
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

% Check for empty epochs and throw error if found
emptyEpochs = find(cellfun(@isempty,{EEG.epoch.event}));
if ~isempty(emptyEpochs)
    error('Found %d empty epochs (indices: %s). This indicates a data processing issue that needs investigation...', ...
          length(emptyEpochs), mat2str(emptyEpochs));
end


end