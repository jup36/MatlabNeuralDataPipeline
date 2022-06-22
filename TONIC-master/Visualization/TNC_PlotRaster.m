function [] = TNC_PlotRaster(figNum,rasterDS)

%% Build up the data array from structure
numTrials = size(rasterDS.trial,2);

for k=1:numTrials

    % create a equal sized trial id vector
    thisTrialTS     = rasterDS.trial(k).ts;
    numStamps       = size(thisTrialTS,1);
    thisTrialIDS    = ones(numStamps,1).*k;

    if k==1
        finalArrayTS = thisTrialTS;            
        finalArrayID = thisTrialIDS;            
    else
        finalArrayTS = [finalArrayTS;thisTrialTS];
        finalArrayID = [finalArrayID;thisTrialIDS];
    end

end

%% Plot routine
figure(figNum);
TNC_CustomRasterPlotter(finalArrayTS,finalArrayID,3,1,[0 0 0],1,1,figNum);

