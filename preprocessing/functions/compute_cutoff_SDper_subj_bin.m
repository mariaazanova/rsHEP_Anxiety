function cutoff_values = compute_cutoff_SDper_subj_bin(EEG_ICA_Corr, initial_multiplier, target_average_ICs, tolerance)
    % Initialize binary search parameters
    max_iterations = 100;
    iter = 0;
    lower_bound = 0;  % Minimum multiplier
    upper_bound = 10; % Maximum multiplier
    multiplier = initial_multiplier;
    
    % Initialize array for cutoff values
    cutoff_values = zeros(length(EEG_ICA_Corr), 1);
    
    % Binary search loop
    while iter < max_iterations
        iter = iter + 1;
        num_removed_ICs = zeros(length(EEG_ICA_Corr), 1);
        
        % Calculate removals with current multiplier
        for i = 1:length(EEG_ICA_Corr)
            subject_mean = mean(EEG_ICA_Corr{i});
            subject_std = std(EEG_ICA_Corr{i});
            cutoff_values(i) = subject_mean + multiplier * subject_std;
            num_removed_ICs(i) = sum(EEG_ICA_Corr{i} > cutoff_values(i));
        end
        
        avg_removed_ICs = mean(num_removed_ICs);
        
        % Check convergence
        if abs(avg_removed_ICs - target_average_ICs) <= tolerance
            fprintf('Converged after %d iterations with multiplier %.3f\n', iter, multiplier);
            break;
        end
        
        % Binary search update
        if avg_removed_ICs > target_average_ICs
            % Too many ICs removed, increase threshold
            lower_bound = multiplier;
            multiplier = (multiplier + upper_bound) / 2;
        else
            % Too few ICs removed, decrease threshold
            upper_bound = multiplier;
            multiplier = (lower_bound + multiplier) / 2;
        end
        
        fprintf('Iteration %d: multiplier=%.3f, avg_ICs=%.2f (target=%.2f)\n', ...
                iter, multiplier, avg_removed_ICs, target_average_ICs);
    end
    
    if iter >= max_iterations
        warning('Maximum iterations reached without convergence');
    end
end