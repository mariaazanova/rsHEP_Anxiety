%% Calculate HEPs
% clc; clear all

% Load the configuration
run('config.m');

%% get files
files = dir(datadir);
files = files([files.isdir] & ~startsWith({files.name}, '.')); % subject folders only

% restrict to the analysis sample (participants remaining after the exclusion criteria)
sample = readtable(sample_file, 'Delimiter', ',', 'TextType', 'string');
sample_ids = unique(strtrim(sample.life_id));
% folder (= raw file) names are LIFE IDs, but two recordings lack the "LI" prefix
folder_ids = string({files.name});
no_prefix = ~startsWith(folder_ids, "LI");
folder_ids(no_prefix) = "LI" + folder_ids(no_prefix);
keep = ismember(folder_ids, sample_ids);
files = files(keep);
life_ids = folder_ids(keep);  

%% documentation
prep_subjects = {};

for hp = 1:length(hp_settings)

    parfor i = 1:length(files)

        try
            % get subjid
            subjid = files(i).name;          
            life_id = char(life_ids(i));     
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
            if ~isfolder(ecgpath)
                mkdir(ecgpath)   
            end

            % load raw EEG and ECG
            [EEG, ECG] = load_life_raw(subjid, raw_dirs, ecg_label, eog_labels, orig_chanlocs);

            %backup of ECG for HEPLAB - only downsampling, rest will be
            %done by heplab
            ECG_heplab = ECG;

            % filter EEG
            % LP at 40 Hz
            EEG = pop_eegfiltnew(EEG, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order], 'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',0,'hicutoff',lp_settings, 'filtorder', [lp_filter_order],'plotfreqz',0);

            % Notch filter
            EEG = pop_eegfiltnew(EEG, 'locutoff',notch_band(1),'hicutoff',notch_band(2),'revfilt',1, 'filtorder', notch_filter_order,'plotfreqz',0);
            ECG = pop_eegfiltnew(ECG, 'locutoff',notch_band(1),'hicutoff',notch_band(2),'revfilt',1, 'filtorder', notch_filter_order,'plotfreqz',0);

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

            % crop to the resting-state block S99 -> S100 (S99 + rest_duration_min if S100 is missing
            % without any start marker the whole recording is kept
            [rest_t0, rest_t1, rest_src] = get_rest_start(EEG, life_id, christian_s99_dirs, rest_duration_min);
            rest_win = [];
            if ~isnan(rest_t0)
                rest_win = [max(rest_t0, 0), min(rest_t1, EEG.pnts / EEG.srate) - 1 / EEG.srate];
                EEG = pop_select(EEG, 'time', rest_win);
                ECG = pop_select(ECG, 'time', rest_win);
                ECG_heplab = pop_select(ECG_heplab, 'time', rest_win);
            end
            EEG.etc.rest_block = struct('start_s', rest_t0, 'end_s', rest_t1, 'source', rest_src, 'kept_min', EEG.pnts / EEG.srate / 60, ...
                                        'win', rest_win);   % [start end] s actually cut; a2 and a9 cut the same window
            ECG.etc.rest_block = EEG.etc.rest_block;

            % ASR for removal of channels (less than correlation of 0.85 or 5 seconds flat)
            EEG = clean_artifacts(EEG, 'ChannelCriterion', chan_crit, 'FlatlineCriterion', flat_line, 'LineNoiseCriterion', 'off', 'BurstCriterion', 'off', 'WindowCriterion', 'off', 'Highpass', 'off');

            % store removed channels in EEG struct
            EEG.reject.removed_channels = all_chan_label(~ismember(all_chan_label , {EEG.chanlocs.labels}));

            % exclude subject if more than max_bad_chan_frac of channels removed
            if numel(EEG.reject.removed_channels) > max_bad_chan_frac * numel(all_chan_label)
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
        
            % R-peaks (2.3.2): HEPLAB, then aligned to the QRS template, R-peaks in unreliable ECG stretches
            % dropped, extrasystoles marked
            HEP_qrs = heplab_slowdetect(ECG_heplab.data, ECG_heplab.srate);
            [n_ext_cers, hr_cers] = get_cers_extrasystoles(cers_file, life_id);
            [HEP_qrs, is_ext, rpeak_info] = process_rpeaks(ECG_heplab.data, HEP_qrs, ECG_heplab.srate, epoch_length, n_ext_cers, hr_cers, rpeak_params);
            parsave([ecgpath, '/', 'extrasystoles.mat'], rpeak_info);

            % Exclude both heartbeats of any pair with an inter-beat interval < min_ibi, and extrasystolic beats
            IBI = diff(HEP_qrs);
            short_ibi = IBI < min_ibi * fs;
            drop_beats = [short_ibi; false] | [false; short_ibi] | is_ext;
            HEP_qrs = HEP_qrs(~drop_beats);

             % save events to EEG and ECG set
            EEG = fun_heplab_save_events(EEG, HEP_qrs, hep_event_marker);
            rs = num2cell(HEP_qrs);
            [EEG.event(end - numel(rs) + 1:end).rpeak_sample] = rs{:};   % R-peak sample in the resting block -> needed for split-half reliability (a10)
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
      
            %% remove all events not at latency 0 (we do not want 2 heartbeats in the same epoch)
            EEG = remove_non_zero_events(EEG);

            % save pre-ICA data - removed baseline from filename
            filename = sprintf('%s_hp_%0.1f.set', subjid, hp_settings(hp));
            pop_saveset(EEG, 'filename', filename, 'filepath', savepath);
            
            ECG_filename = sprintf('%s_hp_%0.1f_ECG.set', subjid, hp_settings(hp));
            pop_saveset(ECG, 'filename', ECG_filename, 'filepath', savepath);
            
            % keep track of processed sujbject IDs 
            prep_subjects{i} = filename;

            %% Some QA Plots
            %pre-clean PSD
            fig = figure('Visible', 'off'); 
            plot_spec(EEG.data(1:length(all_chan_label),:)', 250, 'f_max',70);
            print(fig, [QA_save_path, '/', subjid, '_PSD_pre_ICA.png'], '-dpng', '-r150');
            close(fig);

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