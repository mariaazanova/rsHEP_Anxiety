function cluster_amplitudes = compute_cluster_amplitudes(stat, labelfield, dataFT)
% COMPUTE_CLUSTER_AMPLITUDES Calculate average amplitudes for significant clusters
%
% Inputs:
%   stat       - statistics structure from ft_timelockstatistics
%   labelfield - field name for cluster labels ('posclusterslabelmat' or 'negclusterslabelmat')
%   dataFT     - FieldTrip data structure with individual subject data
%
% Outputs:
%   cluster_amplitudes - structure with amplitude information for each cluster

    % Initialize output structure
    cluster_amplitudes = struct();
    
    % Get unique cluster labels (excluding zeros which indicate non-significant points)
    all_labels = unique(stat.(labelfield));
    all_labels = all_labels(all_labels > 0);
    
    % Process each cluster
    for i = 1:length(all_labels)
        clust_idx = all_labels(i);
        cluster_mask = (stat.(labelfield) == clust_idx);
        
        % Find channels and time points in this cluster
        chan_mask = any(cluster_mask, 2);  % channels in cluster
        time_mask = any(cluster_mask, 1);  % time points in cluster
        
        % Extract data for this cluster
        cluster_data = dataFT.individual(:, chan_mask, time_mask);
        
        % Calculate various measures
        cluster_amplitudes(i).label = clust_idx;
        cluster_amplitudes(i).mean_amplitude = mean(cluster_data(:));
        cluster_amplitudes(i).peak_amplitude = max(abs(mean(mean(cluster_data, 3), 1)));
        cluster_amplitudes(i).std_amplitude = std(cluster_data(:));
        cluster_amplitudes(i).n_channels = sum(chan_mask);
        cluster_amplitudes(i).n_timepoints = sum(time_mask);
        
        % Get time range of cluster
        time_indices = find(time_mask);
        cluster_amplitudes(i).time_start = dataFT.time(time_indices(1));
        cluster_amplitudes(i).time_end = dataFT.time(time_indices(end));
    end
end
