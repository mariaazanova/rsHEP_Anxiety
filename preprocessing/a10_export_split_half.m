%% ROI-average HEP for the first and the second half of the recording 
% Same epochs as a5. 
% Output (derivative_dir): <cond>_split_half_roi_HEP.csv - subject_id, roi_HEP_first, roi_HEP_second, n_first, n_second
% Test runs: set subj_limit before calling; output then goes to <deriv_root>/test/.

run('config.m');

files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); 
out_dir = derivative_dir;
if exist('subj_limit', 'var'), files = files(1:min(subj_limit, end)); out_dir = [deriv_root, 'test/']; end
create_dirs({out_dir});

for ic = 1:length(ic_cutoffs)
    for hp = 1:length(hp_settings)
        for bl = 1:length(baseline)
            if isempty(baseline{bl}), bl_label = 'none'; else, bl_label = sprintf('%d_%d', baseline{bl}); end
            cond_id = sprintf('hp_%0.1f_bl_%s_ic_%0.1f', hp_settings(hp), bl_label, ic_cutoffs(ic));
            res = nan(length(files), 4);                    % first, second, n_first, n_second

            parfor i = 1:length(files)
                try
                    subjid = files(i).name;
                    EEG = pop_loadset('filename', sprintf('%s_hp_%0.1f_ic_%0.1f.set', subjid, hp_settings(hp), ic_cutoffs(ic)), ...
                                      'filepath', fullfile(files(i).folder, subjid, 'clean'));
                    P = pop_loadset('filename', sprintf('%s_hp_%0.1f.set', subjid, hp_settings(hp)), ...
                                    'filepath', fullfile(files(i).folder, subjid, 'prep'), 'loadmode', 'info');
                    if ~isempty(baseline{bl}), EEG = pop_rmbase(EEG, baseline{bl}); end

                    % peak-to-peak rejection and minimum number of epochs, as in a5
                    incl_idx = find(~ismember({EEG.chanlocs.labels}, ignore_channel));
                    [~, ~, bad] = pop_artmwppth_eeglab(EEG, 'Channel', incl_idx, 'Threshold', p2p_threshold, ...
                        'Twindow', epoch_length * 1000, 'Windowsize', p2p_window_ms, 'Windowstep', p2p_step_ms, 'Review', 'off');
                    keep = true(EEG.trials, 1); keep(bad) = false;
                    if sum(keep) < min_epochs, continue, end   % excluded in a5 as well

                    % half of the resting block in which each epoch's R-peak lies
                    onset = nan(EEG.trials, 1);
                    for e = 1:EEG.trials
                        onset(e) = [EEG.event(EEG.epoch(e).event).rpeak_sample];   % only the R-peak event has this field
                    end
                    first = onset <= P.etc.rest_block.kept_min * 60 * EEG.srate / 2;

                    [~, chan_roi] = ismember(electrodes, {EEG.chanlocs.labels});
                    roi = squeeze(mean(mean(EEG.data(chan_roi, EEG.times >= time_range(1) & EEG.times <= time_range(2), :), 1), 2));
                    res(i, :) = [mean(roi(keep & first)), mean(roi(keep & ~first)), sum(keep & first), sum(keep & ~first)];
                catch ME
                    fileID = fopen([save_error, files(i).name, '_', cond_id, '_split_half_error_log.txt'], 'w');
                    fprintf(fileID, '%s\n', getReport(ME));
                    fclose(fileID);
                end
            end

            valid = ~isnan(res(:, 1));
            T = [table({files(valid).name}', 'VariableNames', {'subject_id'}), ...
                 array2table(res(valid, :), 'VariableNames', {'roi_HEP_first', 'roi_HEP_second', 'n_first', 'n_second'})];
            writetable(T, [out_dir, cond_id, '_split_half_roi_HEP.csv']);
        end
    end
end
