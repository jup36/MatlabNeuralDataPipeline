function [binCountsTrialsHz, binEdges] = binTimestampsHz(timeStampC, var)

% Example lick timestamps for 100 trials
%numTrials = 100;
%lick_timestamps_trials = cell(1, numTrials);
%for i = 1:numTrials
%   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
%end

numTrials = length(timeStampC);

% Define bin edges
binEdges = var.minMaxTpre(1):var.binSize:var.minMaxTpost(2); % e.g. from -3s to 8s in 100ms (0.1s) increments

% Initialize storage for binned and z-score normalized lick counts
binCountsTrials = zeros(length(timeStampC), length(binEdges) - 1);

for i = 1:numTrials
    % Bin the lick timestamps for the trial
    binCounts = histc(timeStampC{i}, binEdges);

    % Exclude the last bin
    binCounts = binCounts(1:end-1);
    binCountsTrials(i, :) = binCounts;
end

% convert to Hz
binCountsTrialsHz = binCountsTrials./var.binSize; 

end