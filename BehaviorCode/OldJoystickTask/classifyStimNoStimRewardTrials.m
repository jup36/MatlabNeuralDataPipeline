%This function just loads the BehVariables.mat and classifies the reward
% trials into stim vs no-stim trials, and saves the BehVariables.mat including this index 
function classifyStimNoStimRewardTrials(filePath)
behFile = dir(fullfile(filePath,'BehVariables.mat'));
if length(behFile)>1 || isempty(behFile)    
    error('File could not be found or multiple BehVariables.mat files exist!');
else
    load(behFile.name,'ts') % load behavioral timestamps 
end

rewardStmIdx = logical(cellfun(@(x) sum(abs(x-ts.stmLaser)<=2000), num2cell(ts.reward)));

% just save the index
save(fullfile(filePath,'BehVariables'), 'rewardStmIdx', '-append') % append the behavioral timestamps

end