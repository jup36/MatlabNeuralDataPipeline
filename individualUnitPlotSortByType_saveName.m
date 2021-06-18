function s = individualUnitPlotSortByType_saveName(filePath, dataS, type, psthWin, manualX, saveName)
%This function is to generate an individual unit plot whose trials are
% sorted based on a behavioral variable by the # of fold. 

% put data into a cell - each fold as one cell element
types = unique(type(:,1)); % trial p.Results.trialFolds
typeCell = cell(1,length(types)); 

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 1:length(types)
    typeCell{1,f} = dataS(type(type(:,1)==types(f),2)); % take the trials of the current fold sorted by beh
end
clearvars f 
unitIdFmt = 'Unit#%d'; 
h = spikeRasterGrammSortedFolds( psthWin, manualX, typeCell ); % inputs; psthWin, manualX, typeCell
pbaspect([1 1 1])

s = inputname(3); %strcat(inputname(3),'_',sprintf(unitIdFmt,unit)); 

if isfolder(fullfile(filePath)) % if there's Figure folder already
    print(h,fullfile(filePath,saveName),'-dpdf','-bestfit')
else
    mkdir(fullfile(filePath)) % otherwise, make a folder
    print(h,fullfile(filePath,saveName),'-dpdf','-bestfit')
end

