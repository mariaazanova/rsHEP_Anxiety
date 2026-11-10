%% This script iteratively finds ICA to ECG correlation values that would on average result in 1.5 and 3.5 removed cardiac IC's
% clc; clear all


% Load configuration & settings
run('config.m');

%get files
files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); % subject folders only

for hp = 1:length(hp_settings)
hp_label = sprintf('hp_%0.1f', hp_settings(hp));

%correlation of ICA components with ECG, per subject
EEG_ICA_Corr = cell(1, length(files));
subj_ids = cell(1, length(files));

parfor i=1:length(files)

    try

        %get subjid
        subjid = files(i).name;

        %subject specific folders
        subj_folder = fullfile(files(i).folder, '/',files(i).name, '/');
        eegpath = fullfile(files(i).folder, '/',files(i).name, '/eeg');
        ecgpath = fullfile(files(i).folder, '/',files(i).name, '/ecg');
        ica_path = fullfile(files(i).folder, files(i).name, 'ica');

        %% load ICA data
        ica_id = sprintf('%s_%s_ica', subjid, hp_label);
        EEG_ICA = pop_loadset('filename', [ica_id, '.set'],'filepath',ica_path);
        ECG = pop_loadset('filename', [ica_id, '_ECG.set'], 'filepath', ica_path);

        % %% identify heart components based on correlation with ECG
        %get a template of the ECG to match with IC time courses
        ECG_template =  mean(ECG.data, 3);

        %Select ICA components based on ECG template - SDcrit (1.5) would
        %plot components that have a corrlation with the ECG which is >1.5 standard deviations
        %above the mean correlation across all components. Just for QA 
        [~, cV] = ecg_ica_corr(EEG_ICA, ECG_template, cfa_corr_window, hep_event_marker, 1.5);

        %store ECG ICA correlation in matrix and save it
        EEG_ICA_Corr{i} = cV;
        subj_ids{i} = subjid;

        %List of subjects

        %%if Loop completed and previously a error log was created, delete
        %the error log
        if isfile([save_error, subjid, '_HEP_error_log.txt'])
            delete([save_error, subjid, '_HEP_error_log.txt'])
        end

    catch ME

        fileID = fopen([save_error, subjid, '_', hp_label, '_ICAcutoff_error_log.txt'], 'w');
        fprintf(fileID, '%s\n', getReport(ME));
        fclose(fileID);
    end

end

% keep only subjects that loaded successfully - otherwise threshold computation might fail
ok = ~cellfun(@isempty, subj_ids);
EEG_ICA_Corr = EEG_ICA_Corr(ok);
subj_ids = subj_ids(ok);

save([datadir, 'EEG_ICA_Corr_', hp_label, '.mat'], 'EEG_ICA_Corr')

% Compute cutoff values for each target number of ICs
cutoff_values = struct();
cutoff_sd = struct();
for i = 1:length(target_ICs_removed)
    if target_ICs_removed(i) > 0  
        
        % cutoff label
        cutoff_label = strrep(sprintf('cutoff_value_%.1f', target_ICs_removed(i)), '.', '_'); % 1.5 -> 'cutoff_value_1_5'
        
        [cutoff_values.(cutoff_label), cutoff_sd.(cutoff_label)] = compute_cutoff_SD_per_subj( ...
            EEG_ICA_Corr, cfa_start_sd(i), target_ICs_removed(i), cfa_tolerance);
                
        % Save individual cutoff file
        cutoff_data = cutoff_values.(cutoff_label);
        save([datadir, 'EEG_ICA_', cutoff_label, '_', hp_label, '.mat'],  'cutoff_data');
    end
end

% Save all cutoffs and subject IDs in one file
save([datadir, 'EEG_ICA_cutoffs_', hp_label, '.mat'], 'cutoff_values', 'cutoff_sd');
save([datadir, 'EEG_ICA_cutoff_subj_ids_', hp_label, '.mat'], 'subj_ids');

end