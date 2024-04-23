
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'trI')

%% prepare color
pastel2 = slanCM('Pastel2', 10); 
%plotColorListWithNumbers(pastel2); 

pastel1 = slanCM('Pastel1', 10); 
%plotColorListWithNumbers(pastel1); 

%% control trials classified as low and high load trials 
controlLoI = cellfun(@(x) contains(x, 'lo'), rezCol.controlRsAlignMaxForceId);
controlHiI = cellfun(@(x) contains(x, 'hi'), rezCol.controlRsAlignMaxForceId);

controlForceLo = rezCol.controlRsAlignMaxForce(controlLoI);
controlForceHi = rezCol.controlRsAlignMaxForce(controlHiI); 

%% breakthrough trials classified as left and right trials 
stimFullBtLoI = cellfun(@(x) contains(x, 'lo'), rezCol.stimFullBtForceId);
stimFullBtHiI = cellfun(@(x) contains(x, 'hi'), rezCol.stimFullBtForceId);

stimFullBtForceLo = rezCol.stimFullBtMaxForce(stimFullBtLoI);
stimFullBtForceHi = rezCol.stimFullBtMaxForce(stimFullBtHiI); 

%% stats
[~, stat.maxForceHi.p, ~, stat.maxForceHi.stats] = ttest2(controlForceHi, stimFullBtForceHi);
[~, stat.maxForceLo.p, ~, stat.maxForceLo.stats] = ttest2(controlForceLo, stimFullBtForceLo);




