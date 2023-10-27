function [timeStampOut, poisd, expd] = lickSimPoissonExpBin(timeStampC, minMaxTime, binSize, numReps)
    % Ensure timestamps are within the specified range
    timeStampC = cellfun(@(a) a(a >= minMaxTime(1) & a <= minMaxTime(2)), timeStampC, 'UniformOutput', false);

    % Compute intervals and fit exponential pdf
    intTimeC = cellfun(@(a) diff([minMaxTime(1); a]), timeStampC, 'UniformOutput', false);
    intTime = cell2mat(intTimeC'); 
    expd = fitdist(intTime, 'Exponential');

    % Compute lick counts and fit poisson pdf
    numBins = ceil((minMaxTime(2) - minMaxTime(1)) / binSize);
    binEdges = linspace(minMaxTime(1), minMaxTime(2), numBins + 1);
    counts = cellfun(@(a) histcounts(a, binEdges), timeStampC, 'UniformOutput', false);
    countsMat = cell2mat(counts');
    lambda = mean(countsMat(:)); 
    poisd = makedist('Poisson', 'lambda', lambda);

    % Simulate using the combined model
    timeStampOut = cell(1, numReps);
    for jj = 1:numReps
        simulatedCounts = random(poisd, numBins, 1);
        timeStampOut{jj} = [];
        for bin = 1:numBins
            binStartTime = binEdges(bin);
            timeStampsInBin = binStartTime + cumsum(random(expd, simulatedCounts(bin), 1));
            timeStampOut{jj} = [timeStampOut{jj}; timeStampsInBin(timeStampsInBin < binEdges(bin+1))];
        end
        intervals = diff([minMaxTime(1); timeStampOut{jj}]); 
        intervals = intervals(randperm(length(intervals))); 
        timeStampOut{jj} = minMaxTime(1) +cumsum(intervals); 
    end
end
