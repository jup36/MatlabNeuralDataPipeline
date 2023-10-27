function [lickSimInC, lickSimOutC, lickCntC] = lickSimWrapper(var, rez)
% allocate
lickSimInC = cell(length(var.trSets), 3); % #trSets-by-#epochs, each cell to contain timestamps of # reps
lickSimOutC = cell(length(var.trSets), 3, 2); % #trSets-by-#epochs-by-2, each cell to contain timestamps of # reps
lickCntC = cell(length(var.trSets), 3); % #trSets-by-#epochs, each cell to contain

% divide epochs (pre-, post-, and cue epochs)
for ts = 1:size(lickSimOutC, 1) % increment trial sets (tertiles)
    % pre-cue epoch
    lickSimInC{ts, 1} = cellfun(@(a) a(a < 0), rez.lickTimeC(var.trSets{ts}), 'UniformOutput', false); % pre-cue epoch of the current trial set
    lickCntC{ts, 1} = cellfun(@length, lickSimInC{ts, 1}, 'UniformOutput', false);
    [lickSimOutC{ts, 1, 1}, lickSimOutC{ts, 1, 2}] = lickSimPoissonExpBin(lickSimInC{ts, 1}, var.minMaxTpre, var.simBinSize, var.simRep);

    % cue epoch
    lickSimInC{ts, 2} = cellfun(@(a) a(a >= 0 & a <= var.trDur), rez.lickTimeC(var.trSets{ts}), 'UniformOutput', false); % cue epoch of the current trial set
    lickCntC{ts, 2} = cellfun(@length, lickSimInC{ts, 2}, 'UniformOutput', false);
    [lickSimOutC{ts, 2, 1}, lickSimOutC{ts, 2, 2}] = lickSimPoissonExpBin(lickSimInC{ts, 2}, var.minMaxTcue, var.simBinSize, var.simRep);

    % post-cue epoch
    lickSimInC{ts, 3} = cellfun(@(a) a(a > var.trDur), rez.lickTimeC(var.trSets{ts}), 'UniformOutput', false); % cue epoch of the current trial set
    lickCntC{ts, 3} = cellfun(@length, lickSimInC{ts, 3}, 'UniformOutput', false);
    [lickSimOutC{ts, 3, 1}, lickSimOutC{ts, 3, 2}] = lickSimPoissonExpBin(lickSimInC{ts, 3}, var.minMaxTpost, var.simBinSize, var.simRep);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


end
