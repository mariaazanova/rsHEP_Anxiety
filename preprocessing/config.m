% config.m - LIFE resting-state EEG

%% Paths to tools
addpath /data/p_02865/tools/eeglab2024.1/
addpath /data/p_02865/tools/HEPLAB-master/Functions
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions')) 

orig_chanlocs = '/data/p_02865/tools/eeglab2024.1/plugins/dipfit/standard_BEM/elec/standard_1005.elc';

%% Data paths
raw_dirs = {'/data/p_02865/data/LIFE/raw/EEG/ruhe/', ...
            '/data/p_02865/data/LIFE/raw/EEG/esf_data_2014-08-21/', ...
            '/data/p_02865/data/LIFE/raw/EEG/fehlende_Daten/'};
deriv_root = '/data/p_02865/data/LIFE/derivatives/';
datadir = [deriv_root, 'preprocessing/'];         
save_error = [deriv_root, 'logs/'];
derivative_dir = [deriv_root, 'group/'];
cutoff_dir = [derivative_dir, 'ica_cutoffs/']; 
outdir = [deriv_root, 'output/'];

%ECG and EOG channel labels in the raw data
ecg_label = 'EKG';
eog_labels = {'HEOG', 'VEOG'};

%% Load EEGLAB
eeglab; close;

%in eeglab, option to recompute icaact needs to be "on"; set here, not inside parfor -> causes crashes
pop_editoptions('option_computeica', 1);

%% Pre-processing settings

%HEP event marker
hep_event_marker = '88';

%minimum IBI (sec)
min_ibi = 0.6;

% Length of HEPs in seconds
epoch_length = [-0.2, 0.6]; 

%sampling rate
fs = 250;

%highpass filter cutoff settings
hp_settings = [0.1, 0.5];

%lowpass filter cutoff settings
lp_settings = 40;

%low-pass and notch filters settings
native_fs = 1000;
lp_filter_order = 330;
notch_band = [48, 52];
notch_filter_order = 1650;

%high_pass filter-order settings
hp_filter_order = [8250, 1650];

%high-pass for the separate ICA datasets
ica_hp_cutoff = 1;
ica_hp_filter_order = 826;

%baseline settings
baseline = {[-150, -50], []};

%1.5, 3.5, or no cardiac ICs removed
ic_cutoffs = [1.5, 3.5, 0, NaN];

%clean artifacts settings

%correlation with neighbours
chan_crit = 0.85;

%flat signal duration (seconds)
flat_line = 5;

%exclude participant if more than this % of EEG channels is bad
max_bad_chan_frac = 0.30;

%how many ICs should be removed?
target_ICs_removed = [1.5, 3.5];

%CFA threshold: subject mean + k * SD; k starts here and moves in 0.1 SD steps
%until the average number of removed ICs is within target +/- cfa_tolerance
cfa_start_sd = [3, 1.5]; 
cfa_tolerance = 0.2;
cfa_corr_window = [-0.050, 0.400]; 

%epoch rejection: peak-to-peak amplitude within any moving window
p2p_threshold = 200;
p2p_window_ms = 200;
p2p_step_ms = 100;

%exclude participant from a pipeline if fewer epochs remain after rejection 
min_epochs = 300;

%Ignore channels in bad epoch detection due to influence of eye activity
ignore_channel = {'Fp1', 'Fp2', 'F7', 'F8'};

% region of interest for average HEP calculation
electrodes = {'Cz', 'C3', 'C4', 'Fz', 'F3', 'F4', 'FC1', 'FC2'};
time_range = [200, 500];

% ICA settings
% ICLabel probability thresholds for component rejection
iclabel_threshold_muscle = 0.9;
iclabel_threshold_eye = 0.9;
iclabel_threshold_line_noise = 0.9;
iclabel_threshold_channel_noise = 0.9;
iclabel_threshold_other = 0.9;