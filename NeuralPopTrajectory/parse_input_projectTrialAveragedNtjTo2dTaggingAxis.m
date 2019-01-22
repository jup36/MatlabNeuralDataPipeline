function p = parse_input_projectTrialAveragedNtjTo2dTaggingAxis( filePath, tagPCAFileName, binSpkCntFileName, evtName,  saveNameTag, vargs )

default_binSize = 50; % 50 ms bins
default_kernelSD = 100; % 100 ms sd for the gaussian kernel used for smoothing
default_numbFolds = 3;   % # of folds to divide nTjs into
default_lineWidth  = 1;  % default lineWidth to be used to plot the neural trajectories
default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
default_trjCmap = 'summer'; % the default trajectory colormap
default_PCs = 1:2; % PC directions (axes) to project nTrjs onto

p = inputParser; % create parser object

addRequired(p,'filePath'); % file directory
addRequired(p,'tagPCAFileName'); % fileName for tagging acticity pca result
addRequired(p,'binSpkCntFileName'); % fileName for binSpkCnt
addRequired(p,'saveNameTag'); % saveName used to save the outcomes
addRequired(p,'evtName'); % event name the neural trajectories are aligned to

addParameter(p,'binSize', default_binSize)
addParameter(p,'kernelSD', default_kernelSD)
addParameter(p,'numbFolds', default_numbFolds)
addParameter(p,'lineWidth', default_lineWidth)
addParameter(p,'markerSize', default_markerSize)
addParameter(p,'trjCmap', default_trjCmap)
addParameter(p,'PCs', default_PCs)

parse(p, filePath, tagPCAFileName, binSpkCntFileName, evtName, saveNameTag, vargs{:})
end