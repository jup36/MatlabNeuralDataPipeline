
function [mRez, sRez] = trialGroupMeanSem(tbytPsthC, groupC)
% 'tbytPsthC' contains trial-by-trial data in each cell (it is assumed that data are temporally aligned across trials already).
% 'groupC' contains logical for each grouping, logicals are assumed to be of same lengths as the number of trials.
% 'mRez' returns the group mean for each group logical per row.
% 'sRez' returns the group sem for each group logical per row.

% tbytPsthC = stimOnDff_m1_itp;
% groupC = {[tbytDat.rewardTrI], [tbytDat.punishTrI]};

% sanity check 1: all trial numbers must match!
lenT = length(tbytPsthC);
trialIC = cell2mat(cellfun(@length, groupC, 'UniformOutput', false));
if length(unique([lenT, trialIC]))~=1
    error("The length of trials and length(s) of group logics do not match!")
end

% sanity check 2:
if length(unique(cell2mat(cellfun(@length, tbytPsthC, 'UniformOutput', false))))~=1
    error("Some trials have different lengths, e.g., temporally misaligned!")
end

nGroups = length(groupC);

mRez = [];
sRez = [];

for gg = 1:nGroups
    gDat = cell2mat(tbytPsthC(logical(groupC{gg}), 1));
    [gMean, ~, gSem] = meanstdsem(gDat);
    mRez = [mRez; gMean];
    sRez = [sRez; gSem];
end
end
