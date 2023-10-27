function [timeStampOut, poisd] = lickSimPoisson(timeStampC, minMaxTime, binSize, numReps)
    % Convert timestamps to binned counts
    numBins = ceil((minMaxTime(2) - minMaxTime(1)) / binSize);
    counts = cellfun(@(a) histcounts(a, linspace(minMaxTime(1), minMaxTime(2), numBins + 1)), timeStampC, 'UniformOutput', false);
    countsMat = cell2mat(counts');
    
    % Fit a Poisson distribution
    lambda = mean(countsMat(:)); % average rate across all bins and trials
    poisd = makedist('Poisson', 'lambda', lambda);
    
    % Simulate using the Poisson distribution
    timeStampOut = cell(1, numReps);
    simulatedCounts = cell(1, numReps); 
    for jj = 1:numReps
        simulatedCounts{jj} = random(poisd, numBins, 1);
        timeStampOut{jj} = countsToTimestamps(simulatedCounts{jj}, binSize);
        timeStampOut{jj} = minMaxTime(1) + timeStampOut{jj}(:); 
    end

    %mean(cell2mat(cellfun(@(a) sum(a)/length(a), simulatedCounts, 'un', false)))
    %mean(cell2mat(cellfun(@(a) length(a)./5, timeStampOut, 'un', false)))
    
    
    function timeStamps = countsToTimestamps(counts, binSize)
        timeStamps = [];
        for k = 1:length(counts)
            ts = (k-1)*binSize + binSize * rand(1, counts(k)); % generate timestamps within the bin
            timeStamps = [timeStamps, ts];
        end
    end
end
