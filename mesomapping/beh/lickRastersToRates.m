function [lickRateC] = lickRastersToRates(lickTimestampsC, binEdges)

% Initialize a matrix to store the binned lick counts for all trials
numTrials = length(lickTimestampsC);
binnedLickCounts = zeros(numTrials, length(binEdges) - 1);

% Loop through each trial and calculate the histogram
for trialIdx = 1:numTrials
    % Get the lick timestamps for the current trial
    licks = lickTimestampsC{trialIdx};
    
    % Calculate the histogram
    counts = histcounts(licks, binEdges);
    
    % Store the binned lick counts
    binnedLickCounts(trialIdx, :) = counts;
end

% Convert binned lick counts to lick rates (Hz)
binWidth = binEdges(2) - binEdges(1); % Bin width in seconds
lickRates = binnedLickCounts / binWidth; % Convert counts to rates (licks per second)

end