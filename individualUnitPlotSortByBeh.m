function s = individualUnitPlotSortByBeh(filePath, dataS, sortedBehDat, unit, fold, psthWin, manualX)
%This function is to generate an individual unit plot whose trials are
% sorted based on a behavioral variable by the # of fold. 

% put data into a cell - each fold as one cell element
foldDatCell = cell(1,fold); 
trFolds = floor(length(dataS)/fold); % trial p.Results.trialFolds
% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 0:fold-1
    if f<fold-1
        foldDatCell{1,f+1} = dataS(sortedBehDat(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by beh
    elseif f==fold-1 % take all the remaining trials for the final fold
        foldDatCell{1,f+1} = dataS(sortedBehDat(f*trFolds+1:end,2)); % take the trials of the current fold sorted by beh
    end
end
clearvars f 
unitIdFmt = 'Unit#%d'; 
h = spikeRasterGrammSortedFolds( psthWin, manualX, foldDatCell ); % inputs; psthWin, manualX, foldDatCell
pbaspect([1 1 1])

s = inputname(3); %strcat(inputname(3),'_',sprintf(unitIdFmt,unit)); 

if isfolder(fullfile(filePath,'Figure')) % if there's Figure folder already
    print(h,fullfile(filePath,'Figure',strcat(inputname(3),'_',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
else
    mkdir(fullfile(filePath,'Figure')) % otherwise, make a folder
    print(h,fullfile(filePath,'Figure',strcat(inputname(3),'_',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
end

