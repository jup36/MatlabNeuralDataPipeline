function y = getBinnedSpikeTrain(expt, spLabel, trialIdx)
% y: a sparse column vector representing the concatenated spike trains

sts = cell(numel(trialIdx), 1);
binfun = expt.binfun;
endTrialIndices = [0 cumsum(binfun([expt.trial(trialIdx).duration]))];
nT = endTrialIndices(end); % how many bins total?

for kTrial = trialIdx(:)'
    bst = endTrialIndices(kTrial) + binfun(expt.trial(kTrial).(spLabel));
    sts{kTrial} = bst(:);
end

sts = cell2mat(sts);

y = sparse(sts, 1, 1, nT, 1);
