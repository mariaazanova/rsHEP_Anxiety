%%  copy ICA components and apply ICA cutoff to remove ICs
% clc; clear all

% Load configuration & settings
run('config.m');

%get files
files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); % subject folders only

prep_subjects = {};

for ic = 1:length(ic_cutoffs)

    for hp = 1:length(hp_settings)

        % Load IC cutoffs for this high-pass level
        hp_label = sprintf('hp_%0.1f', hp_settings(hp));
        S = load([datadir, 'EEG_ICA_cutoffs_', hp_label, '.mat']);
        cutoff_values = S.cutoff_values;
        S = load([datadir, 'EEG_ICA_cutoff_subj_ids_', hp_label, '.mat']);
        subj_ids = S.subj_ids;

        parfor i = 1:length(files)

            try

                % Get subjid
                subjid = files(i).name;
                prep_id = sprintf('%s_hp_%0.1f.set', subjid, hp_settings(hp));
                QA_save_id = sprintf('%s_hp_%0.1f_ic_%0.1f', subjid, hp_settings(hp), ic_cutoffs(ic));

                % Subject specific folders
                subj_folder = fullfile(files(i).folder, '/', files(i).name, '/');
                eegpath = fullfile(files(i).folder, '/', files(i).name, '/eeg');
                ecgpath = fullfile(files(i).folder, '/', files(i).name, '/ecg');
                prep_eegpath = fullfile(files(i).folder, '/', files(i).name, '/prep');
                ica_path = fullfile(files(i).folder, files(i).name, 'ica');
                clean_path = fullfile(files(i).folder, files(i).name, 'clean');
                if ~isfolder(clean_path)
                    mkdir(clean_path)
                end
                QA_save_path = fullfile(files(i).folder, '/', files(i).name, '/QA', '/');

                % Load RAW data - only for QA PSD
                EEG_RAW = load_life_raw(subjid, raw_dirs, ecg_label, eog_labels, orig_chanlocs);

                % Load ICA data
                ica_id = sprintf('%s_hp_%0.1f_ica', subjid, hp_settings(hp));
                EEG_ICA = pop_loadset('filename', [ica_id, '.set'], 'filepath', ica_path);

                % Load preprocessed data
                EEG = pop_loadset('filename', prep_id, 'filepath', prep_eegpath);

                % Load ECG 
                ECG = pop_loadset('filename', [ica_id, '_ECG.set'], 'filepath', ica_path);

                % Copy ICA weights to EEG dataset with lower high-pass cutoff
                EEG.etc = EEG_ICA.etc;
                EEG.icawinv = EEG_ICA.icawinv;
                EEG.icasphere = EEG_ICA.icasphere;
                EEG.icaweights = EEG_ICA.icaweights;
                EEG.icachansind = EEG_ICA.icachansind;

                % Again identify heart components based on correlation with ECG
                % Get a template of the ECG to match with IC time courses
                ECG_template = mean(ECG.data, 3);

                % Select ICA components based on ECG template
                [~, cV] = ecg_ica_corr(EEG_ICA, ECG_template, cfa_corr_window, hep_event_marker, 1.5);
               
                %% plot removed ICA components 
                
                % Initialize rej_ic to zeros by default
                rej_ic = zeros(length(cV), 1);
                
                %which is the individual cutoff to take?
                if ~isnan(ic_cutoffs(ic)) && ic_cutoffs(ic) > 0
                    cutoff_label = strrep(sprintf('cutoff_value_%.1f', ic_cutoffs(ic)), '.', '_'); 
                    subj_idx = find(strcmp(subj_ids, subjid));
                    subj_cutoff = cutoff_values.(cutoff_label)(subj_idx);
                    rej_ic = cV > subj_cutoff;
                    plot_ecg_ica_comps(EEG_ICA, rej_ic, cV, QA_save_path, QA_save_id);
                    close all
                end

                 % Heart components
                EEG.etc.ic_remove.heart = find(rej_ic);

                % Remove components from data
                reject = unique([EEG.etc.ic_remove.heart', EEG.etc.ic_remove.muscle', EEG.etc.ic_remove.eye', EEG.etc.ic_remove.LineNoise', EEG.etc.ic_remove.ChannelNoise', EEG.etc.ic_remove.other']);

                %for condition without any ICA component removal, do not remove any components
                if ~isnan(ic_cutoffs(ic))
                    EEGclean = pop_subcomp(EEG, reject);
                else 
                    EEGclean = EEG;
                end
                
                EEGclean = eeg_checkset(EEGclean);

                % Save cleaned dataset - removed baseline from filename
                prep_id_ica = sprintf('%s_hp_%0.1f_ic_%0.1f.set', subjid, hp_settings(hp), ic_cutoffs(ic));
               
                % Saving failed sometimes -> loading ica options in config seems to fix it
                EEGclean = pop_saveset(EEGclean, 'filename', prep_id_ica, 'filepath', clean_path);

                %% Some QA Plots
                try

                %Plot spectogram
                figure; plot_spec(EEGclean.data(1:end,:)', EEGclean.srate, 'f_max',70);
                % pause(2); 
                saveas(gcf, [QA_save_path, QA_save_id, '_PSD_post_ICA.png'],'png');

                %pre-clean PSD
                figure; plot_spec(EEG.data(1:end,:)', EEG.srate, 'f_max',70);
                pause(2); 
                saveas(gcf, [QA_save_path,QA_save_id, '_PSD_pre_ICA.png'],'png');

                %raw PSD
                figure; plot_spec(EEG_RAW.data(1:end,:)', EEG.srate, 'f_max',70);
                pause(2); 
                saveas(gcf, [QA_save_path,QA_save_id, '_PSD_raw.png'],'png');

                %plot removed components
                pop_topoplot(EEG_ICA, 0, reject ,subjid,[] ,0,'electrodes','on','iclabel', 'on', 'iclabel', 'on');
                pause(1); 
                saveas(gcf, [QA_save_path,  QA_save_id, '_removed_components.png'],'png');

                %plot all components
                pop_topoplot(EEG_ICA, 0, 1:size(EEG_ICA.icaweights, 1), subjid, [], 0, 'electrodes', 'on', 'iclabel', 'on');
                pause(5); 

                saveas(gcf, [QA_save_path, QA_save_id, '_all_components.png'],'png');

                %plot artefact components that have high probability to be brain
                %activity
                brain_comps = reject(EEG_ICA.etc.ic_classification.ICLabel.classifications(reject,1)>0.995);
                if any(brain_comps)
                    pop_topoplot(EEG_ICA, 0,brain_comps,subjid,[] ,0,'electrodes','on', 'iclabel', 'on');
                    saveas(gcf, [QA_save_path, QA_save_id, '_potential_brain_components.png'],'png');
                end
                catch ME_plot
                    fileID = fopen([save_error, QA_save_id, '_QAplot_error_log.txt'], 'w');
                    fprintf(fileID, '%s\n', getReport(ME_plot));
                    fclose(fileID);
                end

                close all

                % If loop completed and previously an error log was created, delete the error log
                if isfile([save_error, QA_save_id, '_remICA_error_log.txt'])
                    delete([save_error, QA_save_id, '_remICA_error_log.txt']);
                end
                
            catch ME

                log_id = sprintf('%s_hp_%0.1f_ic_%0.1f', files(i).name, hp_settings(hp), ic_cutoffs(ic));
                fileID = fopen([save_error, log_id, '_remICA_error_log.txt'], 'w');
                fprintf(fileID, '%s\n', getReport(ME));
                fclose(fileID);
            end

        end
    end
end

%% save list of preprocessed subjects
save([datadir, 'prep_subjects_remICA.mat'], 'prep_subjects')