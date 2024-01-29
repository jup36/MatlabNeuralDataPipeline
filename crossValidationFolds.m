function [trainIndices, testIndices] = crossValidationFolds(numTrials, numFolds, trialShuffleLogic)
    % Initialize cell arrays to hold the indices for training and testing sets
    trainIndices = cell(numFolds, 1);
    testIndices = cell(numFolds, 1);

    % Calculate the size of each fold
    foldSize = floor(numTrials / numFolds);

    shuffledTrials = randperm(numTrials);    

    % Assign indices to each fold
    for i = 1:numFolds
        indices = (1 + (i-1)*foldSize) : min(i*foldSize, numTrials);
        
        if trialShuffleLogic
            testIndices{i} = shuffledTrials(indices); 
        else % do not use shuffling
            % Define the testing set for the current fold
            testIndices{i} = indices; 
        end

        % Define the training set as all other indices
        trainIndices{i} = setdiff(1:numTrials, testIndices{i});
    end
end
