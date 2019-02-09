function plotNtrjMultiDimsByTime( filePath, nTrjDimBinTrial, dims, bins, cMap )
%This function takes a 1-folds cell each containing each fold's neural
% trajectories (dimension x timeBins x trials), and plot the PC score of the dimension of interest over time (using gramm pack).

%% plot rank-folded neural population trajectories
x = bins; % x axis

colorMapDims  = cbrewer('div', cMap, length(dims));

figure;
hold on; 
for d = 1:length(dims) % increment dimensions
    dim = dims(d); % the dimension to project onto and plot the projection scores of
    folds_nTrjMat = []; % just to merge across folds by averaging 
    for f = 1:length(nTrjDimBinTrial) % increment folds
        if size(nTrjDimBinTrial{f},3)==1 % there's only one trial or it's already trial-averaged
            fold_nTrjMat = nTrjDimBinTrial{f}(dim,:,:); % get the trial-by-timeBin pcScoreMat of the dim of interest at this fold
        elseif size(nTrjDimBinTrial{f},3)>1 % there are multiple trials, then take the average here
            fold_nTrjMat = permute(squeeze(nTrjDimBinTrial{f}(dim,:,:)),[2 1]); % get the trial-by-timeBin pcScoreMat of the dim of interest at this fold
        end
        folds_nTrjMat = [folds_nTrjMat; fold_nTrjMat]; 
    end    
    [nTrjMatMean,~,nTrjMatSem] = meanstdsem(folds_nTrjMat); % average across trials
    
    if d==2 && length(dims)==2
        boundedline(x,nTrjMatMean,nTrjMatSem,'alpha','transparency',0.1,'cmap',colorMapDims(end,:)); 
    else 
        boundedline(x,nTrjMatMean,nTrjMatSem,'alpha','transparency',0.1,'cmap',colorMapDims(d,:))
    end
end

clearvars f
hold off; 
set(gca,'tickDir','out')
pbaspect([1 1 1])
axis tight
%s=inputname(2); % take the input variable name as a string the name for nTrjDimBinTrial
%figTtlFmt = 'PCAdims#%d'; % figure title format
figTtl = strcat('PCAdims#',num2str(dims)); % format figure title into string
title(figTtl)
legend
print(fullfile(filePath,'Figure',figTtl),'-dpdf');

cd(filePath)
end

