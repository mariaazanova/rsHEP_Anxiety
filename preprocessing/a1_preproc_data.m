%% Calculate HEPs
% clc; clear all

% Load the configuration
run('config.m');

%% get files
files = dir([datadir, 'sub*']);

%% documentation
prep_subjects = {};

for hp = 1:length(hp_settings)

    parfor i = 1:length(files)

        try
            % get subjid
            subjid = files(i).name;
            subj_folder = fullfile(files(i).folder, '/', files(i).name, '/');
            eegpath = fullfile(files(i).folder, '/', files(i).name, '/eeg');
            ecgpath = fullfile(files(i).folder, '/', files(i).name, '/ecg');

            % create output dirs
            savepath = [fullfile(files(i).folder, files(i).name, 'prep')];
            if ~isfolder(savepath)
                mkdir(savepath)
            end

            QA_save_path = [fullfile(files(i).folder, files(i).name, 'QA')];
            if ~isfolder(QA_save_path)
                mkdir(QA_save_path)
            end

            % load data       
            EEG = pop_loadset('filename', [files(i).name, '_task-rest_eeg.set'], 'filepath', eegpath);
            ECG = pop_loadset('filename', [files(i).name, '_task-rest_ecg.set'], 'filepath', ecgpath);
            
            % load channel locations
            EEG = pop_chanedit(EEG, 'lookup', orig_chanlocs);
            
            %Keep only ECG channel
            ecgIdx = find(contains({ECG.chanlocs.labels}, 'ECG'));
            ECG = pop_select(ECG, 'channel', ecgIdx);

            %backup of ECG for HEPLAB - only downsampling, rest will be
            %done by heplab
            ECG_heplab = ECG;

            % filter EEG
            % LP at 40 Hz
            EEG = pop_eegfiltnew(EEG, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order], 'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order],'plotfreqz',0);

            %EEG = pop_eegfilt(EEG, 0, 40, [lp_filter_order], 0, 0, 0, 'fir1', 0);
            %ECG = pop_eegfilt(ECG, 0, 40, [lp_filter_order], 0, 0, 0, 'fir1', 0);

            % Notch filter
            EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',0);

            % downsample
            EEG = pop_resample(EEG, fs);
            ECG = pop_resample(ECG, fs);
            ECG_heplab = pop_resample(ECG_heplab, fs);

            % remove noisy channels
            % keep backup of channels for channel interpolation
            chanlocs = EEG.chanlocs;
            all_chan_label = {chanlocs.labels}

            % HP filter
            EEG = pop_eegfiltnew(EEG, 'locutoff',hp_settings(hp),'hicutoff',0, 'filtorder', [hp_filter_order(hp)],'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',hp_settings(hp),'hicutoff',0, 'filtorder', [hp_filter_order(hp)],'plotfreqz',0);

%           EEG = pop_eegfilt(EEG, hp_settings(hp), 0, hp_filter_order(hp), 0, 0, 0, 'fir1', 0);
%           ECG = pop_eegfilt(ECG, hp_settings(hp), 0, hp_filter_order(hp), 0, 0, 0, 'fir1', 0);

            % ASR for removal of channels (less than correlation of 0.85 or 5 seconds flat)
            EEG = clean_artifacts(EEG, 'ChannelCriterion', chan_crit, 'FlatlineCriterion', flat_line, 'LineNoiseCriterion', 'off', 'BurstCriterion', 'off', 'WindowCriterion', 'off', 'Highpass', 'off');

            % store removed channels in EEG struct
            EEG.reject.removed_channels = all_chan_label(~ismember(all_chan_label , {EEG.chanlocs.labels}));

            % exclude subject if >50% channels removed
            if numel(EEG.reject.removed_channels) >= numel(all_chan_label) / 2
                fileID = fopen([save_error, subjid, '_removed_channel_log.txt'], 'w');
                fprintf(fileID, '%6s\n', subjid);
                fprintf(fileID, '%6s', num2str(numel(EEG.reject.removed_channels)), ' removed channel');
                fclose(fileID);
                continue
            end

            % interpolate
            EEG = pop_interp(EEG, chanlocs, 'spherical');
            eeg_checkset(EEG);

            %average reference
            EEG = pop_reref( EEG, []);
        
            % detect R-peaks from ECG dataset
            HEP_qrs = heplab_slowdetect(ECG_heplab.data, ECG_heplab.srate);

            % %ECG dataset processed in parallel to the EEG
            % ecgIdx = find(contains({ECG.chanlocs.labels}, 'ECG'));
            % ECG = pop_select(ECG, 'channel', ecgIdx);

            % Exclude heartbeats with an inter-beat interval less than 500 ms
            IBI = diff(HEP_qrs); 
            valid_beats = [true; IBI >= min_ibi * fs]; 
            HEP_qrs = HEP_qrs(valid_beats);  

             % save events to EEG and ECG set
            EEG = fun_heplab_save_events(EEG, HEP_qrs, hep_event_marker);
            ECG = fun_heplab_save_events(ECG, HEP_qrs, hep_event_marker);
            eeg_checkset(EEG);
            eeg_checkset(ECG);
            
            % QA: Plot first 30 seconds of ECG with detected R-peaks
            % scale =  abs(min(ECG.data(1,1:30*ECG.srate))-max(ECG.data(1,1:30*ECG.srate)));
            % pop_eegplot(ECG, 1, 1, 1, [], 'winlength', 30, 'spacing', scale);
            % saveas(gcf, [QA_save_path, '/', subjid, '_R_peak_detection_QA.png'], 'png');
            % close;

            % save R-peak onsets
            parsave([ecgpath, '/', 'r_peak_onsets.mat'], HEP_qrs);

            % epoch the data
            EEG = pop_epoch(EEG, {hep_event_marker}, epoch_length);
            ECG = pop_epoch(ECG, {hep_event_marker}, epoch_length);
      
            % save pre-ICA data - removed baseline from filename
            filename = sprintf('%s_hp_%0.1f.set', subjid, hp_settings(hp));
            pop_saveset(EEG, 'filename', filename, 'filepath', savepath);
            
            ECG_filename = sprintf('%s_hp_%0.1f_ECG.set', subjid, hp_settings(hp));
            pop_saveset(ECG, 'filename', ECG_filename, 'filepath', savepath);
            
            % keep track of processed sujbject IDs 
            prep_subjects{i} = filename;

            %% Some QA Plots
            %pre-clean PSD
            figure; plot_spec(EEG.data(1:length(all_chan_label),:)', 250, 'f_max',70);
            saveas(gcf, [QA_save_path, '/', subjid, '_PSD_pre_ICA.png'],'png');

            %% if loop completed and previously an error log was created, delete the error log
            if isfile([save_error, subjid, '_HEP_error_log.txt'])
                delete([save_error, subjid, '_HEP_error_log.txt'])
            end

            % clear EEG ECG

        catch ME
            fileID = fopen([save_error, subjid, '_prep_error_log.txt'], 'w');
            fprintf(fileID, '%6s\n', subjid);
            fprintf(fileID, '%6s\n', getReport(ME));
            fclose(fileID);
        end
    end
end

%% save list of preprocessed subjects
save([datadir, 'prep_subjects.mat'], 'prep_subjects')