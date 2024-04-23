function [trainSetC, testSetC] = balancedResampleTrials(trialLabels, folds)
%This function takes trialLabels and generates train and test sets as cell
% arrays that contain balanced number of each label.

uniqLabels = unique(trialLabels); 

% get random sampled trial numbers that correspond to each label
for i = 1:length(uniqLabels)
    thisLabel = find(trialLabels==uniqLabels(i)); 
    shuffledIndices = randperm(length(thisLabel));
    uniqLabelC{i} = thisLabel(shuffledIndices); 
end

labelNums = cellfun(@length, uniqLabelC); 
minLabelNums = min(labelNums); 
balancedLabelC = cellfun(@(a) a(1:minLabelNums), uniqLabelC, 'UniformOutput', false); 
testSetTrialNum = floor(min(cellfun(@length, uniqLabelC))/folds); 

trainSetC = cell(folds, 1);
testSetC = cell(folds, 1);

% get train and test sets
for i = 1:folds  
    testInd = testSetTrialNum*(i-1)+1:testSetTrialNum*i; 
    testLogic = ismember(1:minLabelNums, testInd)'; 

    testSetC{i, 1} = cell2mat(cellfun(@(a) a(testLogic), balancedLabelC', 'UniformOutput', false)); 
    trainSetC{i, 1} = cell2mat(cellfun(@(a) a(~testLogic), balancedLabelC', 'UniformOutput', false)); 
end

end