function [z_score_binCounts_trials, binEdges] = binTimestampsZscore(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numTrials = length(timeStampC);

        % Define bin edges
        binEdges = var.minMaxTpre(1):var.binSize:var.minMaxTpost(2); % e.g. from -3s to 8s in 100ms (0.1s) increments

        baselineLogic = binEdges < 0;

        % Initialize storage for binned and z-score normalized lick counts
        binCountsTrials = zeros(length(timeStampC), length(binEdges) - 1);

        for i = 1:numTrials
            % Bin the lick timestamps for the trial
            if ~isempty(timeStampC{i})
                binCounts = histc(timeStampC{i}, binEdges);

                % Exclude the last bin
                binCounts = binCounts(1:end-1);
                binCountsTrials(i, :) = binCounts;
            end
        end

        % Calculate mean and standard deviation from the baseline period across all trials (first 50 bins)
        baseline_mean = mean(mean(binCountsTrials(:, baselineLogic)));
        baseline_std = std(binCountsTrials(:, baselineLogic), [], 'all'); % the 'all' option computes the std considering all elements

        % Z-score normalization across all trials
        z_score_binCounts_trials = (binCountsTrials - baseline_mean) / baseline_std;

    end