

function individualUnitPlotSortByType_figSave_each(filePath, dPrm, dataS, type, unit, psthWin, manualX, rank_within_block)
%This function is to generate an individual unit plot whose trials are
% sorted based on a behavioral variable by the # of fold.

% put data into a cell - each fold as one cell element

figSavePath = 'D:\Junchol_Data\JS2p0\collectData\collectFigure\IndividualNeuronPSTH';
wrI = strfind(filePath, 'WR');
fileName = filePath(wrI:wrI+10);

block = num2str(dPrm.maxCoord_2dProj(2));

types = unique(type(:,1)); % trial p.Results.trialFolds
typeCell = cell(1,length(types));

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 1:length(types)
    typeCell{1,f} = dataS(type(type(:,1)==types(f),2)); % take the trials of the current fold sorted by beh
end
clearvars f

unitIdFmt = 'Unit%d';
rankFmt = 'Rank%d';
if dPrm.isStr==0
    saveName = strcat(fileName,'_Ctx','_block', block, '_',sprintf(unitIdFmt, unit), '_', sprintf(rankFmt, rank_within_block));
else
    saveName = strcat(fileName,'_Str','_block', block, '_',sprintf(unitIdFmt, unit), '_', sprintf(rankFmt, rank_within_block));
end

h = spikeRasterGrammSortedFolds( psthWin, manualX, typeCell ); % inputs; psthWin, manualX, typeCell
pbaspect([1 1 1])

print(h,fullfile(figSavePath, strcat('Type',block), saveName), '-dpdf','-bestfit')
close all; 
end
