
load(fullfile('D:\Junchol_Data\JS2p0\collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez')
filePath = 'D:\Junchol_Data\JS2p0\WR40_082019\Matfiles';
load('binSpkCountSTRCTXWR40_082019.mat', 'rStartToPull', 'jkvt', 'spkTimesCell')
S = rStartToPull; 

%%
% individual unit psth aligned to a task event
% sort trials
tq = round([jkvt.pull_torque]./10); % torque
tqT = tq(S.trI)'; tqT(:,2) = 1:length(S.trI); tqT(:,3) = S.trI;  sortByTq = sortrows(tqT,1);
ps = [jkvt.reachP1]; % position
psT = ps(S.trI)'; psT(:,2) = 1:length(S.trI); psT(:,3) = S.trI;  sortByPs = sortrows(psT,1);
tqpsT(:,1) = tqT(:,1)+psT(:,1); tqpsT(:,2) = 1:length(S.trI); tqpsT(:,3) = S.trI; sortByTqPs = sortrows(tqpsT,1); % torque position combo

cellI_stc = 102;
[cellI] = getUnitIdSfromSpkTimesCell(spkTimesCell, S, cellI_stc);
%cellI = 11; % 102 79, 80, 35, 77, 85
thisUnitSpkTimes = S.SpkTimes{cellI};
%individualUnitPlotSortByType(filePath, dPrmRez, thisUnitSpkTimes, sortByTq, cellI, [3e3 2e3], [1e3 1e3]);
%individualUnitPlotSortByType(filePath, dPrmRez, thisUnitSpkTimes, sortByPs, cellI, [3e3 2e3], [1e3 1e3]);
individualUnitPlotSortByType_figSave(filePath, dPrmRez, thisUnitSpkTimes, sortByTqPs, cellI, [3e3 2e3], [1e3 1e3]);






function [unitIdS] = getUnitIdSfromSpkTimesCell(spkTimesCellIn, sIn, unitIdC)
%unit Ids might not match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
% this function identifies the unit from 'S' in 'spkTimesCell' and outputs
% it as 'unitidSTC'

geomI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{4,unitIdC}), sIn.geometry, 'un', 0));
wfI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{6,unitIdC}), sIn.meanWF, 'un', 0));

unitIdS = find(geomI & wfI);

end


function individualUnitPlotSortByType_figSave(filePath, dPrm, dataS, type, unit, psthWin, manualX)
%This function is to generate an individual unit plot whose trials are
% sorted based on a behavioral variable by the # of fold.

% put data into a cell - each fold as one cell element

figSavePath = 'D:\Junchol_Data\JS2p0\collectData\collectFigure\IndividualNeuronPSTH'; 
wrI = strfind(filePath, 'WR');
fileName = filePath(wrI:wrI+10);  

fileI_dPrm = cell2mat(cellfun(@(a) contains(a, fileName), {dPrm.cellId}, 'un', 0)); 
unitI_dPrm = cell2mat(cellfun(@(a) contains(a, ['#',num2str(unit)]), {dPrm.cellId}, 'un', 0)); 
file_unitI = fileI_dPrm & unitI_dPrm; 

block = num2str(dPrm(file_unitI).maxCoord_2dProj(2)); 

types = unique(type(:,1)); % trial p.Results.trialFolds
typeCell = cell(1,length(types));

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 1:length(types)
    typeCell{1,f} = dataS(type(type(:,1)==types(f),2)); % take the trials of the current fold sorted by beh
end
clearvars f

unitIdFmt = 'Unit%d';
saveName = strcat(fileName,'_block', block, '_',sprintf(unitIdFmt, unit));

h = spikeRasterGrammSortedFolds( psthWin, manualX, typeCell ); % inputs; psthWin, manualX, typeCell
pbaspect([1 1 1])

print(h,fullfile(figSavePath, strcat('Type',block), saveName), '-dpdf','-bestfit')

end
