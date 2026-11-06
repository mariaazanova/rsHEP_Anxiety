%Selection of ICA components based on correlation with ECG Template
function [EEG, cV, rejV] = ecg_ica_corr(EEG, ECG_template, time_window, R_marker, SDcrit)
% Inputs:
%   EEG          - EEGLAB EEG structure with ICA components
%   ECG_template - ECG template for correlation (optional)
%   time_window  - Time window around R peaks for epoching [start, end] in seconds
%   R_marker     - Marker for R peaks in EEG event structure
%   SDcrit       - Standard deviation criterion for cardiac component selection
%
% Outputs:
%   EEG   - EEGLAB EEG structure after processing
%   cV    - Correlation values between components and ECG template
%   rejV  - Logical vector indicating components marked for rejection



EEGtemp = EEG;

%in eeglab, option to recompute icaact needs to be "on"
pop_editoptions('option_computeica', 1)
EEGtemp = eeg_checkset(EEGtemp)

%overwrite the EEG data with the ICA time data for easier epoching / data
%handling
EEGtemp.data = EEGtemp.icaact;

%if data is not epched yet, do it here
if ndims(EEGtemp.data) <=2

    %exchange EEG data with ICA data for easier epoching
    if size(EEGtemp.data,1) ~= size(EEG.data,1);
        add_chans =  size(EEG.data,1) -size(EEGtemp.data,1)
        EEGtemp.data = [EEGtemp.data; zeros(add_chans, length(EEGtemp.data))];
    end

    %Epoch ICA data
    EEGtemp = pop_epoch(EEGtemp, {R_marker} , time_window);

    EEGtemp.data(end-add_chans:end,:,:) = [];

    %Epoch EEG according to ECG marker
    EEG = pop_epoch(EEG, {R_marker} , time_window); 
end

% average ICA time courses
avgIC = mean(EEGtemp.data,3);

%if no ECG template was given, try to create one
if isempty(ECG_template)
    try
        ECG_template = mean(EEG.data(strcmp({EEG.chanlocs.labels},'ECG'),:,:),3);
    catch
    end
end

%Correlate with the ECG template
cV  = abs(corr(avgIC',ECG_template')); 

%obtain threshold
corthreshV = mean(cV)+ SDcrit * std(cV); 

%select comps which are over threshold
rejV = cV > corthreshV; 

%component indices that could be CFA
comps = find(rejV);

if sum(comps) == 0
    disp('No ECG components identified')
else
    % 
    % %check if oscillatory peak is present
    % for k=1:sum(rejV)
    % 
    %     comp = comps(k);
    %     alpha_rel_pow =  Compute_RP(EEGtemp.data(comp,:), 500,[1,45], [7 15]);
    %     %                 disp([num2str(comp), ' has relative alpha pow of  ', num2str(alpha_rel_pow)])
    %     if alpha_rel_pow > 0.5  % if peak greater than 50% - keep the component
    %         rejV(comps(k)) = 0;
    %     else
    %     end
    % 
    % end
    % 
end



end
