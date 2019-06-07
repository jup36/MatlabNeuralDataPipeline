function rewardTrialsWithOutStim(filePath)
%This script is just to sort out reward trials with or without stim trials

cd(filePath)
load(fullfile(filePath,'BehVariables.mat'),'ts')

rewardTrialRange = arrayfun(@(x) x-2000:x, ts.reward', 'UniformOutput', false ); 
latestStimPerReward = arrayfun(@(x) ts.stmLaser(find(ts.stmLaser<x, 1, 'last')), ts.reward', 'UniformOutput', false); 

stmRwdIdx = zeros(length(ts.reward),1); 
for t=1:length(ts.reward)
    if ~isempty(latestStimPerReward{t})
        stmRwdIdx(t) = ismember(latestStimPerReward{t},rewardTrialRange{t}); 
    end
end

ts.rewardNoStim = ts.reward(~stmRwdIdx); 

save(fullfile(filePath,'BehVariables.mat'),'ts','-append')

end

