%% Obtain average HEPs per subject and calculate time/space average HEP

% Load configuration & settings
run('config.m');

%get files
files = dir([datadir, 'sub*']);

for ic = 1:length(ic_cutoffs)
    for hp = 1:length(hp_settings)

        for bl = 1:length(baseline)
        
            % Initialize local variables with NaNs for parfor
            avg_HEP_local = nan(length(files), 1);
            HEP_all_local = cell(length(files), 1);
            ECG_all_local = cell(length(files), 1);

            parfor i = 1:length(files)
                try

                    % Get subjid
                    subjid = files(i).name;
                    QA_save_id = sprintf('%s_hp_%0.1f_bl_%s_ic_%0.1f', subjid, hp_settings(hp), mat2str(baseline{bl}), ic_cutoffs(ic));
                    
                    % Load ICA-cleaned data (without baseline in filename)
                    prep_id_ica = sprintf('%s_hp_%0.1f_ic_%0.1f.set', subjid, hp_settings(hp), ic_cutoffs(ic));
                    prep_id_ica_ECG = sprintf('%s_hp_%0.1f_ECG.set', subjid, hp_settings(hp));

                    % Subject specific folders
                    subj_folder = fullfile(files(i).folder, '/', files(i).name, '/');
                    eegpath = fullfile(files(i).folder, '/', files(i).name, '/eeg');
                    ecgpath = fullfile(files(i).folder, '/', files(i).name, '/ecg');
                    prep_eegpath = fullfile(files(i).folder, '/', files(i).name, '/prep');
                    HEP_path = fullfile(files(i).folder, '/', files(i).name, '/HEP/');
                    QA_save_path = fullfile(files(i).folder, '/', files(i).name, '/QA/');

                    if ~isfolder(HEP_path)
                        mkdir(HEP_path)
                    end

                    % load cleaned dataset
                    EEG = pop_loadset('filename', prep_id_ica, 'filepath', prep_eegpath);
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

                    [~, ~, bad_epoch_indices, ptp_amplitudes] = pop_artmwppth_eeglab(EEG, 'Channel', incl_idx, 'Threshold', 200, 'Twindow', [-200 800], 'Windowsize', 200, 'Windowstep', 100, 'Review', 'off');

                    % remove the epochs
                    if ~isempty(bad_epoch_indices)
                        EEG = pop_rejepoch(EEG, bad_epoch_indices, 0);
                        EEG.reject.ptp.epochs = bad_epoch_indices;
                        EEG.reject.ptp.amplitudes = ptp_amplitudes;

                        % Store rejection statistics
                        rejected_epochs_info(i).subject = subjid;
                        rejected_epochs_info(i).n_total = n_epochs_pre;
                        rejected_epochs_info(i).n_rejected = length(bad_epoch_indices);
                        rejected_epochs_info(i).percent_rejected = 100 * length(bad_epoch_indices) / n_epochs_pre;
                
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
                    avg_HEP_local(i) = mean(mean(meanHEP(chan_roi, time_roi)));

                catch
                    % Failed entries remain NaN
                end
            end

            % Remove invalid entries using NaN check
            valid_idx = ~isnan(avg_HEP_local);
            avg_HEP_local = avg_HEP_local(valid_idx);

            % Get and save valid subject IDs
            all_subject_ids = {files.name};
            valid_subject_ids = all_subject_ids(valid_idx);
            valid_subjects_file = sprintf('hp_%0.1f_bl_%s_ic_%0.1f_valid_subjects.txt', hp_settings(hp), mat2str(baseline{bl}), ic_cutoffs(ic));
            full_path = [derivative_dir valid_subjects_file];
            writecell(valid_subject_ids', full_path, 'Delimiter', 'tab');

            % Convert tabs to newlines using system command
            system(['sed -i "s/\t/\n/g" "' full_path '"']);

            % Use same valid indices for HEP and ECG data
            HEP_all_local = HEP_all_local(valid_idx);
            HEP_all = cat(3, HEP_all_local{:});
            ECG_all_local = ECG_all_local(valid_idx);
            ECG_all = cat(3, ECG_all_local{:});

            % Save with baseline info in filename
            all_hep_id = sprintf('hp_%0.1f_bl_%s_ic_%0.1f_all_HEP.mat', hp_settings(hp), mat2str(baseline{bl}), ic_cutoffs(ic));      
            all_ecg_id = sprintf('hp_%0.1f_bl_%s_ic_%0.1f_all_ECG.mat', hp_settings(hp), mat2str(baseline{bl}), ic_cutoffs(ic));      

            save([derivative_dir, all_hep_id], 'HEP_all');
            save([derivative_dir, all_ecg_id], 'ECG_all');

        end
    end
end

% QA plots for rejected epochs
figure;
hist([rejected_epochs_info.percent_rejected], 20);
title('Distribution of Rejected Epochs Percentage');
xlabel('Percentage of Epochs Rejected');
ylabel('Number of Subjects');
saveas(gcf, [derivative_dir 'QA_rejected_epochs_distribution.png']);
close;
