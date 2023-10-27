function [timeStampOut, poisd, expd] = lickSimPoissonExp(timeStampC, minMaxTime, binSize, numReps)
% Ensure timestamps are within the specified range
timeStampC = cellfun(@(a) a(a >= minMaxTime(1) & a <= minMaxTime(2)), timeStampC, 'UniformOutput', false);

% Compute intervals and fit exponential pdf
intTimeC = cellfun(@(a) diff([minMaxTime(1); a]), timeStampC, 'UniformOutput', false);
intTime = cell2mat(intTimeC');
expd = fitdist(intTime, 'Exponential');

% Compute lick counts and fit poisson pdf
numBins = ceil((minMaxTime(2) - minMaxTime(1)) / binSize);
counts = cellfun(@(a) histcounts(a, linspace(minMaxTime(1), minMaxTime(2), numBins + 1)), timeStampC, 'UniformOutput', false);
countsMat = cell2mat(counts');
lambda = mean(countsMat(:));
poisd = makedist('Poisson', 'lambda', lambda);

% Simulate using the combined model
timeStampOut = cell(1, numReps);
for jj = 1:numReps
    simulatedCounts = random(poisd, numBins, 1);
    timeStampOut{jj} = countsToTimestamps(simulatedCounts, minMaxTime, expd);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function timeStamps = countsToTimestamps(simulCounts, minMaxTime, expd)
        timeStamps = [];
        intervals = random(expd, sum(simulCounts), 1);
        intSum = cumsum(intervals);

        if isempty(intSum)
            timeStamps = [];
        else
            if intSum(end) > diff(minMaxTime)
                intervals = intervals(intSum <= diff(minMaxTime));

                if isempty(intervals)
                    intervals = 0;
                else
                    for i = 1:100
                        intervals = [intervals; random(expd, 1)];
                        intSum = cumsum(intervals);
                        intervals = intervals(intSum <= diff(minMaxTime));
                    end
                    if length(intervals) > length(simulCounts)
                        intervals = intervals(1:length(simulCounts)); 
                    end
                end
            end
            timeStamps = minMaxTime(1) + intSum;
            timeStamps = timeStamps(timeStamps < minMaxTime(2));
        end
    end
end
