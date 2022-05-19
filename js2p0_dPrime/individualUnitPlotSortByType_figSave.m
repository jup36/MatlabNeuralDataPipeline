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
