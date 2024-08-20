function zScoreNormalized = zscoreNormalizationUnitTimeTrial(spikeCounts, numBaseBins) 

[numUnit, numBin, numTr] = size(spikeCounts); 

% Assuming your tensor is named 'spikeCounts'
% Extract the baseline period 
baselinePeriod = spikeCounts(:, 1:numBaseBins, :);

% Reshape the baselinePeriod to collapse the time and trial dimensions
baselinePeriodReshaped = reshape(baselinePeriod, [numUnit, numBaseBins * numTr]);

% Calculate the mean and standard deviation across all time bins and trials for each neuron
meanBaseline = mean(baselinePeriodReshaped, 2); % Mean across the combined time and trial dimension
stdBaseline = std(baselinePeriodReshaped, 0, 2); % Standard deviation across the combined time and trial dimension

% Replicate the mean and standard deviation to match the size of the original tensor
meanBaseline = repmat(meanBaseline, [1, numBin, numTr]);
stdBaseline = repmat(stdBaseline, [1, numBin, numTr]);

% Perform z-score normalization
zScoreNormalized = (spikeCounts - meanBaseline) ./ stdBaseline;

% Handling cases where stdBaseline is zero to avoid NaN values
zScoreNormalized(stdBaseline == 0) = NaN;


end