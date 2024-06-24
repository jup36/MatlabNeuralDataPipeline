function dividedTrials = divideTrialsIntoBlocks(totalTrials, numSets, numTrials, excludeStart, excludeEnd)
    % Function to divide trials into specified number of sets
    % excluding some trials from the start and end.

    if nargin < 3
        excludeStart = 0;
    end

    if nargin < 4
        excludeEnd = 0;
    end

    % Remove specified number of trials from start and end
    trimmedTrials = (excludeStart + 1):(totalTrials - excludeEnd);


    % Divide trials into sets
    dividedTrials = cell(1, numSets);
    for i = 1:numSets
        startIdx = (i-1) * numTrials + 1;
        endIdx = min(i * numTrials, totalTrials); % Ensure we don't exceed the trial count
        trialIdx = startIdx:endIdx; 

        dividedTrials{i} = trialIdx(ismember(trialIdx, trimmedTrials));
    end
end
