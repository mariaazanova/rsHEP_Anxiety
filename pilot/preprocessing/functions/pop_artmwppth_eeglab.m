function [EEG, com, bad_epoch_indices, ptp_amplitudes] = pop_artmwppth_eeglab(EEG, varargin)
% POP_ARTMWPPTH_EEGLAB - Artifact rejection based on moving window peak-to-peak threshold
%
% Usage:
%   >> [EEG, com, bad_epoch_indices, ptp_amplitudes] = pop_artmwppth_eeglab(EEG, 'key', value, ...);
%
% Inputs:
%   EEG         - EEGLAB EEG structure
%
% Optional inputs (key-value pairs):
%   'Twindow'   - Time window for analysis in milliseconds [start_ms end_ms] (default: [EEG.xmin*1000 EEG.xmax*1000])
%   'Threshold' - Peak-to-peak amplitude threshold in microvolts (default: 100)
%   'Windowsize'- Size of the moving window in milliseconds (default: 200)
%   'Windowstep'- Step size for the moving window in milliseconds (default: 100)
%   'Channel'   - Channels to include in the analysis (default: all channels)
%   'LowPass'   - Low-pass filter cutoff frequency in Hz (default: -1, no filtering)
%   'Flag'      - Flag to mark rejected epochs (default: 1)
%   'Review'    - Option to review data ('on' or 'off', default: 'off')
%
% Outputs:
%   EEG                - EEGLAB EEG structure after artifact rejection
%   com                - Command string for EEGLAB history
%   bad_epoch_indices  - Indices of epochs marked as bad
%   ptp_amplitudes     - Peak-to-peak amplitudes of the rejected epochs
%
% Example:
%   [EEG, com, bad_epochs, amplitudes] = pop_artmwppth_eeglab(EEG, 'Threshold', 75, 'Channel', [1:32]);


com = '';
bad_epoch_indices = [];  % Initialize the array to store bad epoch indices
ptp_amplitudes = [];     % Initialize the array to store the peak-to-peak amplitudes

if nargin < 1
    help pop_artmwppth_eeglab;
    return;
end

% Parsing inputs
p = inputParser;
p.FunctionName  = mfilename;
p.CaseSensitive = false;
p.addRequired('EEG');
p.addParameter('Twindow', [EEG.xmin*1000 EEG.xmax*1000], @isnumeric);
p.addParameter('Threshold', 100, @isnumeric);
p.addParameter('Windowsize', 200, @isnumeric);
p.addParameter('Windowstep', 100, @isnumeric);
p.addParameter('Channel', 1:EEG.nbchan, @isnumeric);
p.addParameter('LowPass', -1, @isnumeric);
p.addParameter('Flag', 1, @isnumeric);
p.addParameter('Review', 'off', @ischar);

p.parse(EEG, varargin{:});

testwindow =  p.Results.Twindow;
ampth      =  p.Results.Threshold;
winms      =  p.Results.Windowsize;
stepms     =  p.Results.Windowstep;
chanArray  =  p.Results.Channel;
lpval      =  p.Results.LowPass;
flag       =  p.Results.Flag;

fs       = EEG.srate;
winpnts  = floor(winms*fs/1000);
stepnts  = floor(stepms*fs/1000);

% Error checking
if stepnts < 1
    error('Step value should be equal to or greater than the sampling period (1/fs msec).');
end

% Filtering if required
if lpval > 1
    EEG_lowfilt = pop_eegfiltnew(EEG, [], lpval, 8260, 0, [], 0);
else
    EEG_lowfilt = EEG;
end

% Artifact rejection process
max_ptp = zeros(1, EEG.trials);  % Initialize array to store max ptp amplitude for each epoch
for ch = 1:length(chanArray)
    chan = chanArray(ch);
    for epoch = 1:EEG.trials
        for j = 1:stepnts:(EEG.pnts-winpnts+1)
            window_data = EEG_lowfilt.data(chan, j:j+winpnts-1, epoch);
            vs = abs(max(window_data) - min(window_data));
            if vs > ampth
                max_ptp(epoch) = max(max_ptp(epoch), vs);  % Update max ptp for this epoch
                EEG.reject.rejmanual(epoch) = true;
                EEG.reject.rejmanualE(chan, epoch) = true;
            end
        end
    end
end

% Collect the indices and max ptp amplitudes of the rejected epochs
bad_epoch_indices = find(EEG.reject.rejmanual);
ptp_amplitudes = max_ptp(bad_epoch_indices)';

% If Review option is enabled
if strcmpi(p.Results.Review, 'on')
    pop_eegplot(EEG, 1, 1, 1);
end

EEG = eeg_checkset(EEG);

% Command history
com = sprintf('%s = pop_artmwppth_eeglab(%s, ''Twindow'', [%d %d], ''Threshold'', %d, ''Windowsize'', %d, ''Windowstep'', %d, ''Channel'', %s, ''LowPass'', %d, ''Flag'', %d, ''Review'', ''%s'');', ...
    inputname(1), inputname(1), testwindow(1), testwindow(2), ampth, winms, stepms, mat2str(chanArray), lpval, flag, p.Results.Review);

% Display completion
fprintf('Artifact rejection completed.\n');
fprintf('Bad epoch indices: %s\n', mat2str(bad_epoch_indices));  % Display the indices of the bad epochs
disp('Peak-to-peak amplitudes of rejected epochs:');
disp(ptp_amplitudes);
end
