%% Obtain average HEPs aross subjects
% clc; clear all

% Load configuration & settings
run('config.m');

%% load matrices contaning all per subject average HEP / ECG
hep_files = dir([derivative_dir, '*all_HEP.mat*'])
ecg_files = dir([derivative_dir, '*all_ECG.mat*'])

% load single dataset just for layout and time vector
% any cleaned dataset of this run works as the template
tmpl = dir(fullfile(datadir, '*', 'clean', '*_hp_0.5_ic_3.5.set'));
EEG = pop_loadset('filename', tmpl(1).name, 'filepath', tmpl(1).folder);

%% region of interest
electrodes = {'Cz', 'C3', 'C4', 'Fz', 'F3', 'F4', 'FC1', 'FC2'}
time_range = [200, 500]

%% get average across space and time hep
[~, chan_roi] = ismember(electrodes, {EEG.chanlocs.labels});
time_roi = EEG.times >= time_range(1) & EEG.times <= time_range(2);

%% store all ROI average HEPs to plot in single plot
grand_all_HEP = [];
conditions = struct();

for i = 1:length(hep_files)

    try
        
        %load data
        load([derivative_dir, hep_files(i).name])
        load([derivative_dir, ecg_files(i).name])

        %calculate average across trials
        avg_HEP = mean(HEP_all,3);
        avg_ECG = mean(ECG_all,3);

        % Add info about which condition this represents
        [~,fname,~] = fileparts(hep_files(i).name);
        parts = split(fname,'_');
        hp = str2double(parts{2});
        bl = strjoin(parts(4:5),'_');
        ic = str2double(parts(7));
        
        conditions(i).hp = hp;
        conditions(i).baseline = bl;
        conditions(i).ic = ic;

        %To plot all HEPs
        grand_all_HEP(i,:) = mean(avg_HEP(chan_roi,:));
        
    catch ME
        warning('Error processing file %s: %s', hep_files(i).name, ME.message);
    end
end

% save avg frontal ROI HEP
fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600]); 
hold on
for plt = 1:size(grand_all_HEP,1)
    plot(EEG.times, grand_all_HEP(plt,:), 'DisplayName', ...
        sprintf('hp=%.1f bl=%s ic=%.1f', ...
        conditions(plt).hp, conditions(plt).baseline, conditions(plt).ic))
end
legend('show', 'Location', 'northeastoutside','Interpreter', 'none');
xlabel('time (ms)')
ylabel('mikrovolt')
title('Grand Average HEP over frontal electrodes')
print(fig, [outdir, 'frontal_all_HEP.png'], '-dpng', '-r150');
close(fig)
