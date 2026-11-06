function [cutoff_values, multiplier] = compute_cutoff_SD_per_subj(EEG_ICA_Corr, initial_multiplier, target_average_ICs, tolerance)
% Subject mean + k * subject SD, one k for the sample; k moves in 0.1 SD steps until the
% sample-average number of removed ICs is within target +/- tolerance 
% Errors if no k is found 

    % Initialize the multiplier
    multiplier = initial_multiplier;

    % Initialize an array to store the cutoff values for each subject
    cutoff_values = zeros(length(EEG_ICA_Corr), 1);

    % Iterate until the average number of removed ICs is within the desired range
    for iter = 1:200   
        num_removed_ICs = zeros(length(EEG_ICA_Corr), 1);  % To store the number of removed ICs for the current threshold

        % Calculate the cutoff and number of ICs exceeding the cutoff for each subject
        for i = 1:length(EEG_ICA_Corr)
            % Calculate the mean and standard deviation for the current subject
            subject_mean = mean(EEG_ICA_Corr{i});
            subject_std = std(EEG_ICA_Corr{i});

            % Set the initial cutoff for this subject
            cutoff_values(i) = subject_mean + multiplier * subject_std;

            % Calculate the number of ICs exceeding the cutoff
            num_removed_ICs(i) = sum(EEG_ICA_Corr{i} > cutoff_values(i));
        end

        % Calculate the average number of removed ICs across subjects
        avg_removed_ICs = mean(num_removed_ICs);

        % Check if the average number of removed ICs falls within the target range
        if abs(avg_removed_ICs - target_average_ICs) <= tolerance
            return
        else
            % Adjust the multiplier if the average number of removed ICs is outside the range
            if avg_removed_ICs < target_average_ICs - tolerance
                multiplier = multiplier - 0.1;  % Decrease the multiplier to include more ICs
            elseif avg_removed_ICs > target_average_ICs + tolerance
                multiplier = multiplier + 0.1;  % Increase the multiplier to include fewer ICs
            end
        end
    end
    error('No multiplier within target %.1f +/- %.1f ICs after 200 steps (last k = %.1f, avg = %.2f ICs)', ...
          target_average_ICs, tolerance, multiplier, avg_removed_ICs);
end
