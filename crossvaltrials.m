function [ testTrialLogic ] = crossvaltrials( totalTrial, k )
%This function returns trainTrialLogic to select train trials with replacement 
% among the totalTrials for k-fold crossvalidation. 

% Sanity check
if floor(totalTrial/k)<1
    error('Check the inputs for crossvaltrials!')
end

seed = (1:totalTrial)';      % total trials to be sampled
testTrialLogic = cell(k,1);  % train trial index
used = zeros(totalTrial,1);  % to mark the test trials already used (this is to ensure testing all trials)

% get test trial sets for k-fold cross validation
for i = 1:k  % increment k folds
    tmpTestTrials      = zeros(totalTrial,1); % initialize train trials
    if i < k    
        testTrials         = datasample(seed(~used), floor(totalTrial/k), 'Replace', false); % these will be test trials
        used(testTrials,1) = true;         % mark the used trials for test
        tmpTestTrials(testTrials,1) = 1;   % get test trials for this set
        testTrialLogic{i,1}= logical(tmpTestTrials);
        testTrialLogic{i,2}= sort(testTrials);
    elseif i == k % at the final iteration
        testTrials         = seed(~used); % include all the unused trials
        tmpTestTrials(testTrials,1) = 1;   % get test trials for this set
        testTrialLogic{i,1}= logical(tmpTestTrials);
        testTrialLogic{i,2}= sort(testTrials);
    end
end

end

