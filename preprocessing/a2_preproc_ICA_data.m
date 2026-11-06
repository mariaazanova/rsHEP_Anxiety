%% Calculate HEPs
% clc; clear all

% Load configuration & settings
run('config.m');

%get files
files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); % subject folders only

%cell array for preprocessed subject IDs
prep_subjects = {}

for hp = 1:length(hp_settings)

    parfor i = 1:length(files)

        try

            %get subjid
            subjid = files(i).name;

            %preprocessed subjid - take one as a reference for removed channels
            prep_id = sprintf('%s_hp_%0.1f.set', subjid, hp_settings(hp));
            
            %subject specific folders
            subj_folder = fullfile(files(i).folder, '/',files(i).name, '/');
            eegpath = fullfile(files(i).folder, '/',files(i).name, '/eeg');
            ecgpath = fullfile(files(i).folder, '/',files(i).name, '/ecg');
            prep_eegpath = fullfile(files(i).folder, '/',files(i).name, '/prep');

            % output dir
            savepath = fullfile(files(i).folder, files(i).name, 'ica');
            if ~isfolder(savepath)
                mkdir(savepath)
            end

            %% load RAW data - will be filtered at 1 Hz and used for ICA
            % ECG is loaded raw again as the one from a1_preproc_data.m is already high-pass filtered and epoched
            [EEG_ICA, ECG] = load_life_raw(subjid, raw_dirs, ecg_label, eog_labels, orig_chanlocs);

            %% load preprocessed dataset from subject
            EEG = pop_loadset('filename', prep_id,'filepath', prep_eegpath);

            %% filter EEG
            %Low-Pass ICA filter
            EEG_ICA = pop_eegfiltnew(EEG_ICA, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order], 'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order],'plotfreqz',0);

            %% Notch filter, identical to a1_preproc_data.m (before downsampling)
            EEG_ICA = pop_eegfiltnew(EEG_ICA, 'locutoff',notch_band(1),'hicutoff',notch_band(2),'revfilt',1, 'filtorder', notch_filter_order,'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',notch_band(1),'hicutoff',notch_band(2),'revfilt',1, 'filtorder', notch_filter_order,'plotfreqz',0);

            %% downsample
            EEG_ICA = pop_resample(EEG_ICA, fs);
            ECG = pop_resample(ECG, fs);

            %% Remove noisy channels
            chanlocs = EEG.chanlocs;
            all_chan_label = {chanlocs.labels};

            %remove channels identified in previousfrom ICA data
            if isfield(EEG.reject, 'removed_channels')
                EEG_ICA = pop_select(EEG_ICA, 'nochannel', EEG.reject.removed_channels); %remove ECG channel
                EEG_ICA.reject.removed_channels = EEG.reject.removed_channels;
                rmChan = numel(EEG.reject.removed_channels);
            else
                EEG_ICA.reject.removed_channels = [];
                rmChan = 0;
            end

            %% High pass filter
            EEG_ICA = pop_eegfiltnew(EEG_ICA, 'locutoff',ica_hp_cutoff,'hicutoff',0, 'filtorder', ica_hp_filter_order,'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',ica_hp_cutoff,'hicutoff',0, 'filtorder', ica_hp_filter_order,'plotfreqz',0);

            %% Interpolate
            EEG_ICA = pop_interp(EEG_ICA, chanlocs, 'spherical');
            eeg_checkset(EEG_ICA);

            %% Average reference, identical to the target datasets in a1_preproc_data.m
            EEG_ICA = pop_reref(EEG_ICA, []);

            %% Add R-peaks detected and IBI-cleaned in a1_preproc_data.m
            HEP_qrs = parload([ecgpath, '/', 'r_peak_onsets.mat']);
        
            EEG_ICA = fun_heplab_save_events(EEG_ICA, HEP_qrs, hep_event_marker);
            ECG = fun_heplab_save_events(ECG, HEP_qrs, hep_event_marker);
            eeg_checkset(EEG_ICA);

            %% Epoch the data
            EEG_ICA = pop_epoch(EEG_ICA, {hep_event_marker} , epoch_length);
            ECG = pop_epoch(ECG, {hep_event_marker} , epoch_length);

            %% remove all events not at latency 0 (we do not want 2 heartbeats in the same epoch)
            EEG_ICA = remove_non_zero_events(EEG_ICA);

            % %% remove some epochs for testing
            % EEG_ICA = pop_selectevent( EEG_ICA, 'event',[1:20] ,'deleteevents','off','deleteepochs','on','invertepochs','off');

            %% run extended infomax ICA
            %get number of non-interpolated channels
            % ind_channel =  length(EEG_ICA.chanlocs)-rmChan;

            % rank of the data
            dataRank = sum(eig(cov(double(EEG_ICA.data(:,:)'))) > 1E-7)

            EEG_ICA = pop_runica(EEG_ICA, 'icatype', 'runica', 'pca', dataRank , 'options', {'extended' 1});


            %% run IClabel for classification of components
            EEG_ICA = pop_iclabel(EEG_ICA, 'default');

            %Add above threshold component idx to data
            EEG_ICA.etc.ic_remove.muscle = find(EEG_ICA.etc.ic_classification.ICLabel.classifications(:,2) >= iclabel_threshold_muscle);
            EEG_ICA.etc.ic_remove.eye = find(EEG_ICA.etc.ic_classification.ICLabel.classifications(:,3) >= iclabel_threshold_eye);
            EEG_ICA.etc.ic_remove.LineNoise = find(EEG_ICA.etc.ic_classification.ICLabel.classifications(:,5) >= iclabel_threshold_line_noise);
            EEG_ICA.etc.ic_remove.ChannelNoise = find(EEG_ICA.etc.ic_classification.ICLabel.classifications(:,6) >= iclabel_threshold_channel_noise);
            EEG_ICA.etc.ic_remove.other = find(EEG_ICA.etc.ic_classification.ICLabel.classifications(:,7) >= iclabel_threshold_other);

            %% Save ICA EEG and corresponding ECG data
            ica_id = sprintf('%s_hp_%0.1f_ica', subjid, hp_settings(hp));
            pop_saveset(EEG_ICA, 'filename', [ica_id, '.set'], 'filepath', savepath);
            pop_saveset(ECG, 'filename', [ica_id, '_ECG.set'], 'filepath', savepath);

            %list of all preprocessed subjects
            prep_subjects{i} = subjid;

            %if Loop completed and previously a error log was created, delete
            %the error log
            if isfile([save_error, subjid, '_error_log_ICA_comp_select.txt'])
                delete([save_error, subjid, '_error_log_ICA_comp_select.txt'])
            end

            close all

            %%if Loop completed and previously a error log was created, delete
            %the error log
            if isfile([save_error, subjid, '_HEP_error_log.txt'])
                delete([save_error, subjid, '_HEP_error_log.txt'])
            end

        catch ME

            fileID = fopen([save_error, subjid, '_ICA_error_log.txt'],'w');
            fprintf(fileID,'%6s\n',subjid);
            fprintf(fileID,'%6s\n',getReport(ME));
            fclose(fileID);
        end

    end
end

%% save list of preprocessed subjects
save([datadir, 'prep_subjects_ICA.mat'], 'prep_subjects')