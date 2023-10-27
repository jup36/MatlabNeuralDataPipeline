function dividedTrials = divideTrials(totalTrials, numSets, excludeStart, excludeEnd)
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

    % Calculate set lengths
    numTrials = length(trimmedTrials);
    setLength = ceil(numTrials / numSets);

    % Divide trials into sets
    dividedTrials = cell(1, numSets);
    for i = 1:numSets
        startIdx = (i-1) * setLength + 1;
        endIdx = min(i * setLength, numTrials); % Ensure we don't exceed the trial count
        
        dividedTrials{i} = trimmedTrials(startIdx:endIdx);
    end
end
