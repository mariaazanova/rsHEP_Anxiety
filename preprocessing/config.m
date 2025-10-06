% config.m

%% Paths to tools
addpath /data/p_02865/tools/eeglab2024.1/
addpath /data/p_02865/tools/HEPLAB-master/Functions
addpath /data/p_02865/code/preprocessing/functions

%% Data paths
datadir = '/data/p_02865/data/ds003838/';
save_error = '/data/p_02865/data/ds003838/Logfiles/';
derivative_dir = [datadir, 'derivatives/'];
outdir = [datadir, 'output/'];

orig_chanlocs = '/data/p_02865/tools/eeglab2024.1/plugins/dipfit/standard_BEM/elec/standard_1005.elc';

%% Load EEGLAB
eeglab; close;

%in eeglab, option to recompute icaact needs to be "on"

%% Pre-processing settings

%HEP event marker
hep_event_marker = '88';

%minimum IBI (sec)
min_ibi = 0.5;

% Length of HEPs in seconds
epoch_length = [-0.2, 0.8]

%sampling rate
fs = 250

%highpass filter cutoff settings
hp_settings = [0.1, 0.5];

%lowpass filter cutoff settings
lp_settings = 40;

%lowpass filter-order settings
lp_filter_order = 86; %changed to 86 so its even and filtnew can be used

%high_pass filter-order settings
hp_filter_order = [8250, 1650];

%baseline settings
baseline = {[-150, -50], []};

%1.5, 3.5, or no cardiac ICs removed
ic_cutoffs = [1.5, 3.5, 0, NaN];

%clean artifacts settings

%correlation with neighbours
chan_crit = 0.85;

%flat signal duration (seconds)
flat_line = 5;

%how many ICs should be removed?
target_ICs_removed = [1.5, 3.5, 0];


%Ignore channels in bad epoch detection due to influence of eye activity
ignore_channel = {'Fp1', 'Fp2', 'F7', 'F8'};

% region of interest for average HEP calculation
electrodes = {'Cz', 'C3', 'C4', 'Fz', 'F3', 'F4', 'FC1', 'FC2'};
time_range = [200, 500];