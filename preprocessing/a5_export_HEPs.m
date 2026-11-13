%% Obtain average HEPs per subject and calculate time/space average HEP

% Load configuration & settings
run('config.m');

%get files
files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); % subject folders only

create_dirs({derivative_dir, [qa_dir, 'epoch_rejection/'], [qa_dir, 'valid_subjects/']});

for ic = 1:length(ic_cutoffs)
    for hp = 1:length(hp_settings)

        for bl = 1:length(baseline)
        
            % Initialize local variables with NaNs for parfor
            avg_HEP_local = nan(length(files), 1);
            HEP_all_local = cell(length(files), 1);
            ECG_all_local = cell(length(files), 1);
            n_total_local = nan(length(files), 1);     
            n_rejected_local = nan(length(files), 1);  

            % construct name
            if isempty(baseline{bl}), bl_label = 'none'; else, bl_label = sprintf('%d_%d', baseline{bl}); end
            cond_id = sprintf('hp_%0.1f_bl_%s_ic_%0.1f', hp_settings(hp), bl_label, ic_cutoffs(ic));

            parfor i = 1:length(files)
                try

                    % Get subjid
                    subjid = files(i).name;
                    QA_save_id = [subjid, '_', cond_id];
                    
                    % Load ICA-cleaned data (without baseline in filename)
                    prep_id_ica = sprintf('%s_hp_%0.1f_ic_%0.1f.set', subjid, hp_settings(hp), ic_cutoffs(ic));
                    prep_id_ica_ECG = sprintf('%s_hp_%0.1f_ECG.set', subjid, hp_settings(hp));

                    % Subject specific folders
                    subj_folder = fullfile(files(i).folder, '/', files(i).name, '/');
                    eegpath = fullfile(files(i).folder, '/', files(i).name, '/eeg');
                    ecgpath = fullfile(files(i).folder, '/', files(i).name, '/ecg');
                    prep_eegpath = fullfile(files(i).folder, '/', files(i).name, '/prep');
                    clean_path = fullfile(files(i).folder, files(i).name, 'clean');
                    HEP_path = fullfile(files(i).folder, '/', files(i).name, '/HEP/');
                    QA_save_path = fullfile(files(i).folder, '/', files(i).name, '/QA/');

                    if ~isfolder(HEP_path)
                        mkdir(HEP_path)
                    end

                    % load cleaned dataset
                    EEG = pop_loadset('filename', prep_id_ica, 'filepath', clean_path);
                    ECG = pop_loadset('filename', prep_id_ica_ECG, 'filepath', prep_eegpath);

                    % apply baseline correction
                    if ~isempty(baseline{bl})
                        EEG = pop_rmbase(EEG, baseline{bl});
                        ECG = pop_rmbase(ECG, baseline{bl});
                    end

                    %% peak-to-peak epoch removal
                    incl_idx = find(~ismember({EEG.chanlocs.labels}, ignore_channel));

                    % identify bad epochs exceeding some peak-to-peak
                    % amplitude thresholds
                    n_epochs_pre = EEG.trials;

                    [~, ~, bad_epoch_indices, ptp_amplitudes] = pop_artmwppth_eeglab(EEG, 'Channel', incl_idx, 'Threshold', p2p_threshold, 'Twindow', epoch_length * 1000, 'Windowsize', p2p_window_ms, 'Windowstep', p2p_step_ms, 'Review', 'off');

                    % remove the epochs
                    if ~isempty(bad_epoch_indices)
                        EEG = pop_rejepoch(EEG, bad_epoch_indices, 0);
                        ECG = pop_rejepoch(ECG, bad_epoch_indices, 0);
                        EEG.reject.ptp.epochs = bad_epoch_indices;
                        EEG.reject.ptp.amplitudes = ptp_amplitudes;
                    end
                    
                    % Store rejection statistics 
                    n_total_local(i) = n_epochs_pre;
                    n_rejected_local(i) = length(bad_epoch_indices);

                    % minimum number of remaining epochs
                    if EEG.trials < min_epochs
                        error('Only %d epochs left after rejection (minimum %d)', EEG.trials, min_epochs);
                    end

                    %% remove all events not at latency 0 (we do not want 2 heartbeats in the same epoch)
                    EEG = remove_non_zero_events(EEG);

                    %% calculate average across trials
                    meanHEP = mean(EEG.data, 3);
                    meanECG = mean(ECG.data, 3);

                    % Store data from all subjects
                    HEP_all_local{i} = meanHEP;
                    ECG_all_local{i} = meanECG;

                    %% get average across space and time hep
                    [~, chan_roi] = ismember(electrodes, {EEG.chanlocs.labels});
                    time_roi = EEG.times >= time_range(1) & EEG.times <= time_range(2);

                    % Store average HEP (now using NaN as invalid marker)
                    roi_HEP = mean(mean(meanHEP(chan_roi, time_roi)));
                    avg_HEP_local(i) = roi_HEP;

                    % subject wise average HEP/ECG for this condition 
                    subj_HEP = struct('meanHEP', meanHEP, 'meanECG', meanECG, 'times', EEG.times, ...
                        'chanlabels', {{EEG.chanlocs.labels}}, 'roi_HEP', roi_HEP, ...
                        'n_epochs_total', n_epochs_pre, 'n_epochs_rejected', length(bad_epoch_indices), ...
                        'n_epochs_kept', size(EEG.data, 3));
                    parsave([HEP_path, subjid, '_', cond_id, '_avgHEP.mat'], subj_HEP);

                catch ME
                    % Failed entries remain NaN -> the error per pipeline
                    fileID = fopen([save_error, QA_save_id, '_HEP_error_log.txt'], 'w');
                    fprintf(fileID, '%s\n', getReport(ME));
                    fclose(fileID);
                end
            end

            % Remove invalid entries using NaN check
            valid_idx = ~isnan(avg_HEP_local);

            % epoch rejection per subject for this condition
            has_data = ~isnan(n_total_local);
            rej_table = table({files(has_data).name}', n_total_local(has_data), n_rejected_local(has_data), ...
                n_total_local(has_data) - n_rejected_local(has_data), valid_idx(has_data), ...
                'VariableNames', {'subject_id', 'n_total', 'n_rejected', 'n_kept', 'included'});
            writetable(rej_table, [qa_dir, 'epoch_rejection/', cond_id, '_epoch_rejection.csv']);

            % QA: distribution of rejected epochs for this condition
            fig = figure('Visible', 'off');
            histogram(100 * rej_table.n_rejected ./ rej_table.n_total, 20);
            title(strrep(cond_id, '_', ' ')); xlabel('% epochs rejected'); ylabel('Number of subjects');
            print(fig, [qa_dir, 'epoch_rejection/', cond_id, '_QA_rejected_epochs.png'], '-dpng', '-r150');
            close(fig);
            avg_HEP_local = avg_HEP_local(valid_idx);

            % Get and save valid subject IDs
            all_subject_ids = {files.name};
            valid_subject_ids = all_subject_ids(valid_idx);
            valid_subjects_file = [cond_id, '_valid_subjects.txt'];
            full_path = [qa_dir, 'valid_subjects/', valid_subjects_file];
            writecell(valid_subject_ids', full_path, 'Delimiter', 'tab');

            % Convert tabs to newlines using system command
            system(['sed -i "s/\t/\n/g" "' full_path '"']);

            % ROI-average HEP per valid subject (same subject order as HEP_all)
            roi_table = table(valid_subject_ids', avg_HEP_local, n_total_local(valid_idx), ...
                n_rejected_local(valid_idx), n_total_local(valid_idx) - n_rejected_local(valid_idx), ...
                'VariableNames', {'subject_id', 'roi_HEP', 'n_epochs_total', 'n_epochs_rejected', 'n_epochs_kept'});
            writetable(roi_table, [derivative_dir, cond_id, '_roi_HEP.csv']);

            % Use same valid indices for HEP and ECG data
            HEP_all_local = HEP_all_local(valid_idx);
            HEP_all = cat(3, HEP_all_local{:});
            ECG_all_local = ECG_all_local(valid_idx);
            ECG_all = cat(3, ECG_all_local{:});

            % Save with baseline info in filename
            all_hep_id = [cond_id, '_all_HEP.mat'];
            all_ecg_id = [cond_id, '_all_ECG.mat'];

            subject_ids = valid_subject_ids';
            save([derivative_dir, all_hep_id], 'HEP_all', 'subject_ids');
            save([derivative_dir, all_ecg_id], 'ECG_all', 'subject_ids');

        end
    end
end
