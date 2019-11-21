function s = individualUnitPlotSortByType(filePath, dataS, type, unit, psthWin, manualX)
%This function is to generate an individual unit plot whose trials are
% sorted based on a behavioral variable by the # of fold. 

% put data into a cell - each fold as one cell element
types = unique(type(:,1)); % trial p.Results.trialFolds
typeCell = cell(1,length(types)); 

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 1:length(types)
    typeCell{1,f} = dataS(type(:,1)==types(f)); % take the trials of the current fold sorted by beh
end
clearvars f 
unitIdFmt = 'Unit#%d'; 
h = spikeRasterGrammSortedFolds( psthWin, manualX, typeCell ); % inputs; psthWin, manualX, typeCell
pbaspect([1 1 1])

s = inputname(3); %strcat(inputname(3),'_',sprintf(unitIdFmt,unit)); 

if isfolder(fullfile(filePath,'Figure')) % if there's Figure folder already
    print(h,fullfile(filePath,'Figure',strcat(inputname(3),'_',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
else
    mkdir(fullfile(filePath,'Figure')) % otherwise, make a folder
    print(h,fullfile(filePath,'Figure',strcat(inputname(3),'_',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
end

