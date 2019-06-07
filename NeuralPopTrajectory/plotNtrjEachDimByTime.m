function plotNtrjEachDimByTime( filePath, nTrjDimBinTrial, dim, bins, cMap )
%This function takes a 1-folds cell each containing each fold's neural
% trajectories (dimension x timeBins x trials), and plot the PC score of the dimension of interest over time (using gramm pack).

%% plot rank-folded neural population trajectories
x = bins; % x axis
% get trajectories of the current dimension for each trial fold to be trial-by-timeBins
%y = cellfun(@(x) permute(squeeze(x(dim,:,:)),[2 1]), nTrjDimBinTrial, 'UniformOutput', false)'; % trajectories permuted to be trial-by-timeBins
% get class lables
nTrjMat = []; 
c = [];

if ischar(cMap) % extract colors per fold
    colorMapFolds = TNC_CreateRBColormapJP(length(nTrjDimBinTrial),cMap);
elseif size(cMap,1)==length(nTrjDimBinTrial) && size(cMap,2)==3 % the input might be a color map per fold
    colorMapFolds = cMap; 
end

figure;
hold on; 
for f = 1:length(nTrjDimBinTrial) % increment folds
    if size(nTrjDimBinTrial{f},3)==1 % there's only one trial or it's already trial-averaged
        fold_nTrjMat = nTrjDimBinTrial{f}(dim,:,:); % get the trial-by-timeBin pcScoreMat of the dim of interest at this fold     
    elseif size(nTrjDimBinTrial{f},3)>1 % there are multiple trials, then take the average here 
        fold_nTrjMat = permute(squeeze(nTrjDimBinTrial{f}(dim,:,:)),[2 1]); % get the trial-by-timeBin pcScoreMat of the dim of interest at this fold         
    end
    [fold_nTrjMatMean,~,fold_nTrjMatSem] = meanstdsem(fold_nTrjMat); 
    boundedline(x,fold_nTrjMatMean,fold_nTrjMatSem,'alpha','transparency',0.1,'cmap',colorMapFolds(f,:))
    nTrjMat = [nTrjMat; fold_nTrjMat]; % accumulate the nTrjMat
    % get class labels
    cFold = sprintf('%01d',zeros(1,size(fold_nTrjMat,1))+f);
    c = [c;cFold'];
end
clearvars f
hold off; 
set(gca,'tickDir','out')
pbaspect([1 1 1])
axis tight
s=inputname(2); % take the input variable name as a string the name for nTrjDimBinTrial
figTtlFmt = 'PCAdim#%d'; % figure title format
figTtl = sprintf(figTtlFmt,dim); % format figure title into string
title(strcat(s,figTtl)); 

cd(fullfile(filePath,'Figure'))
print(strcat(s,figTtl),'-dpdf');


%c = cellstr(c);
%nTrjMatCell = mat2cell(nTrjMat,ones(size(nTrjMat,1),1),size(nTrjMat,2)); % convert the nTrjMat to a cell array

% clear g
% g(1,1)=gramm('x',x,'y',nTrjMatCell,'color',c);
% g(1,1).stat_summary();
% g(1,1).set_title('stat_summary()');
% 
% g.set_title('');
% figure('Position',[100 100 200 200]);
% g.draw();
% clearvars dim y c

cd(filePath)
end

