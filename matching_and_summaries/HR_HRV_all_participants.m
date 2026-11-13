
%% Calculate HR and HRV (RMSSD) for all LIFE participants

clear; clc;

% Paths and parameters
base_dir = '/data/p_02865/data/LIFE/derivatives/preprocessing';
output_file = fullfile(fileparts(base_dir), 'HR_HRV_all_participants.csv');
fs = 250;  % ECG sampling frequency (Hz)

% Get participant folders
folders = dir(base_dir);
folders = folders([folders.isdir]);
folders = folders(~startsWith({folders.name}, '.'));

nSub = numel(folders);

% Initialize results
subid = strings(nSub, 1);
HR_bpm = nan(nSub, 1);
RMSSD_ms = nan(nSub, 1);

%% Loop through participants
for i = 1:nSub

    subid(i) = string(folders(i).name);

    % Path to R-peak events
    file_path = fullfile(base_dir, subid(i), ...
        'ecg', 'r_peak_onsets.mat');

    if ~isfile(file_path)
        warning('Missing R-peaks: %s', subid(i));
        continue;
    end

    try
        % Load R-peak sample indices
        S = load(file_path);
        r_peaks = double(S.x(:));

        % Require at least 3 peaks for RMSSD
        if numel(r_peaks) < 3 || ...
                any(~isfinite(r_peaks)) || ...
                any(diff(r_peaks) <= 0)
            warning('Invalid R-peaks: %s', subid(i));
            continue;
        end

        % Interbeat intervals (milliseconds)
        RR_ms = diff(r_peaks) / fs * 1000;

        % Mean heart rate (beats per minute)
        HR_bpm(i) = 60000 / mean(RR_ms);

        % HRV: Root Mean Square of Successive Differences
        RMSSD_ms(i) = sqrt(mean(diff(RR_ms).^2));

    catch ME
        warning('Error for %s: %s', subid(i), ME.message);
    end

end

%% Save results
results = table(subid, HR_bpm, RMSSD_ms);

writetable(results, output_file);

fprintf('Saved HR/HRV for %d participants to:\n%s\n', ...
    sum(~isnan(RMSSD_ms)), output_file);
