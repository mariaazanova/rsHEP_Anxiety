%wrapper script to run HEP analysis 
% clc; clear all;

% Run each script sequentially
scripts = {
    'a1_preproc_data.m',
    'a2_preproc_ICA_data.m',
    'a3_find_ICA_cutoff.m',
    'a4_remove_ICA_components.m',
    'a5_export_HEPs.m',
    'a6_grand_avg_HEPs.m'
};


for i = 1:length(scripts)
    run(scripts{i});
end
