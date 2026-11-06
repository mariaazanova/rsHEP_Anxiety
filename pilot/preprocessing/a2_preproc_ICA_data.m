%% Calculate HEPs
% clc; clear all

% Load configuration & settings
run('config.m');

%get files
files = dir([datadir, 'sub*']);

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
            savepath = [fullfile(files(i).folder, files(i).name, 'prep')];

            %% load RAW data - will be filtered at 1Hz and used for ICA
            EEG_ICA = pop_loadset('filename', [files(i).name, '_task-rest_eeg.set'],'filepath',eegpath);

            %ECG is epoched at the end of prep_a1 - hence applying a
            %different HP filter might not be a good idea. Hence lets
            %load the raw ECG again for the ICA processing
            ECG = pop_loadset('filename', [files(i).name, '_task-rest_ecg.set'],'filepath',ecgpath);

            %Keep only ECG channel
            ecgIdx = find(contains({ECG.chanlocs.labels}, 'ECG'));
            ECG = pop_select(ECG, 'channel', ecgIdx);

            %% load  preprocessed dataset from subject
            EEG = pop_loadset('filename', prep_id,'filepath', prep_eegpath);

            %% filter EEG
            %Low-Pass ICA filter
            EEG_ICA = pop_eegfiltnew(EEG_ICA, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order], 'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order],'plotfreqz',0);

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

            %% Interpolate
            EEG_ICA = pop_interp(EEG_ICA, chanlocs, 'spherical');
            eeg_checkset(EEG_ICA);

            %% High pass filter
            EEG_ICA = pop_eegfiltnew(EEG_ICA, 'locutoff',1,'hicutoff',0, 'filtorder', [hp_filter_order(hp)],'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',1,'hicutoff',0, 'filtorder', [hp_filter_order(hp)],'plotfreqz',0);
            
            %% Notch filter
            EEG_ICA = pop_eegfiltnew(EEG_ICA, 'locutoff',48, 'hicutoff',52,'revfilt',1,'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',0);

            %% Add ECG markers previously obtained and saved
            % HEP_qrs = parload([ecgpath, '/', 'r_peak_onsets.mat']);

            % detect R-peaks in ECG
            HEP_qrs = heplab_slowdetect(ECG.data, ECG.srate);

            EEG_ICA = fun_heplab_save_events(EEG_ICA, HEP_qrs, hep_event_marker);
            ECG = fun_heplab_save_events(ECG, HEP_qrs, hep_event_marker);
            eeg_checkset(EEG_ICA);

            %% Epoch the data
            EEG_ICA = pop_epoch(EEG_ICA, {hep_event_marker} , epoch_length);
            ECG = pop_epoch(ECG, {hep_event_marker} , epoch_length);

            %% remove all events not at latency 0 (we do not want 2 heartbeats in the same epoch)
            EEG = remove_non_zero_events(EEG);

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
            pop_saveset(EEG_ICA, 'filename', [subjid, '_ica'], 'filepath', savepath);
            pop_saveset(ECG, 'filename', [subjid, '_ica_ECG'], 'filepath', savepath);

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