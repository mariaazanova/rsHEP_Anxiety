%% This script iteratively finds ICA to ECG correlation values that would on average result in 1.5 and 3.5 removed cardiac IC's
% clc; clear all


% Load configuration & settings
run('config.m');

%get files
files = dir([datadir, 'sub*']);

%matrix for correlation of ICA components with ECG
EEG_ICA_Corr = [];

%struct for successfully processed subjects
prep_subjects = {};

counter = 1;
for i=1:length(files)

    try

        %get subjid
        subjid = files(i).name;

        %subject specific folders
        subj_folder = fullfile(files(i).folder, '/',files(i).name, '/');
        eegpath = fullfile(files(i).folder, '/',files(i).name, '/eeg');
        ecgpath = fullfile(files(i).folder, '/',files(i).name, '/ecg');
        prep_eegpath = fullfile(files(i).folder, '/',files(i).name, '/prep');

        %% load ICA data
        EEG_ICA = pop_loadset('filename', [subjid, '_ica.set'],'filepath',prep_eegpath);
        ECG = pop_loadset('filename', [subjid, '_ica_ECG.set'], 'filepath', prep_eegpath);

        % %% identify heart components based on correlation with ECG
        %get a template of the ECG to match with IC time courses
        ECG_template =  mean(ECG.data, 3);

        %Select ICA components based on ECG template - SDcrit (1.5) would
        %plot components that have a corrlation with the ECG which is >1.5 standard deviations
        %above the mean correlation across all components. Just for QA 
        [EEG, cV, rejV] = ecg_ica_corr(EEG_ICA,ECG_template, [-0.050 0.400], hep_event_marker, 1.5);

        %store ECG ICA correlation in matrix and save it
        EEG_ICA_Corr{counter} = cV;
        subj_ids{counter} = subjid;
        counter = counter+1;

        %List of subjects

        %%if Loop completed and previously a error log was created, delete
        %the error log
        if isfile([save_error, subjid, '_HEP_error_log.txt'])
            delete([save_error, subjid, '_HEP_error_log.txt'])
        end

    catch ME

        fileID = fopen([save_error, subjid, '_remICA_error_log.txt'],'w');
        fprintf(fileID,'%6s\n',subjid);
        fprintf(fileID,'%6s\n',getReport(ME));
        fclose(fileID);
    end

end

save([datadir, 'EEG_ICA_Corr.mat'], 'EEG_ICA_Corr')

% Compute cutoff values for each target number of ICs
cutoff_values = struct();
for i = 1:length(target_ICs_removed)
    if target_ICs_removed(i) > 0  
        
        % cutoff label
        cutoff_label = sprintf('cutoff_value_%d_%d', ...
            floor(target_ICs_removed(i)), ...
            round(10 * mod(target_ICs_removed(i), 1)));
        
        cutoff_values.(cutoff_label) = compute_cutoff_SDper_subj_bin(EEG_ICA_Corr, ...
            target_ICs_removed(i), target_ICs_removed(i), 0.2);
                
        % Save individual cutoff file
        cutoff_data = cutoff_values.(cutoff_label);
        save([datadir, 'EEG_ICA_', cutoff_label, '.mat'],  'cutoff_data');
    end
end

% Save all cutoffs and subject IDs in one file
save([datadir, 'EEG_ICA_cutoffs.mat'], 'cutoff_values');
save([datadir, 'EEG_ICA_cutoff_subj_ids.mat'], 'subj_ids');