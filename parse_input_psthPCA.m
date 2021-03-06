function p = parse_input_psthPCA( filePath, fileName, varName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
%parse input, and extract name-value pairs for the main function
% 'pcaPSTH.m'

default_binSize  = 200; % 200 ms bins
default_stepSize = default_binSize; % by default, use the stepSize same as the binSize
default_binSizeZ = 50;  % binSize for SpkCountMatZ
default_PCs      = 3;   % use top 3 PCs
default_expVarCut = 80; % to include top PCs whose summed explained variance greater than 80% of total variance
default_FRcut = 1; % exclude units with mean FR lower than the FRcut from PCA
default_nanTrialCut = .1; % to exclude NaN trials, if the proportion of them relative to total trials is smaller than this proportion
default_nanUnitCut = .1;  % to exclude units with NaN trials, if the proportion of them relative to total units is smaller than this proportion
default_useAllUnits = true;  % logical to include all units in the sortedBinSpkCell prepared for DCA
default_PCsLogic = false;    % logical to choose to include only the top PCs units
default_expVarLogic = false; % logical to choose to include units with PCs of which cumulative explained variance surpassing the set expVarCut (e.g. 80 %)
default_cmap = 'cb';    % default colormap
default_cAxis = [-3 3]; % default colorAxis
default_rmvNaNunits = false; % exclude units with NaN trials from PCA and DCA preprocessing
default_rmvNaNtrials = true; % remove the trials in which one or more units have NaN trials
default_imagescXlim = [-2000 2000]; % default xlim to be used for imagesc plot 1ms 

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'fileName');
addRequired(p,'varName');

addParameter(p,'binSize', default_binSize)
addParameter(p,'stepSize',default_stepSize)
addParameter(p,'binSizeZ', default_binSizeZ)
addParameter(p,'PCs', default_PCs)
addParameter(p,'expVarCut', default_expVarCut)
addParameter(p,'FRcut', default_FRcut)
addParameter(p,'nanTrialCut', default_nanTrialCut)
addParameter(p,'nanUnitCut', default_nanUnitCut)
addParameter(p,'useAllUnits',default_useAllUnits)
addParameter(p,'PCsLogic',default_PCsLogic)
addParameter(p,'expVarLogic',default_expVarLogic)
addParameter(p,'cmap',default_cmap)
addParameter(p,'cAxis',default_cAxis)
addParameter(p,'rmvNaNunits',default_rmvNaNunits)
addParameter(p,'rmvNaNtrials',default_rmvNaNtrials)
addParameter(p,'imagescXlim',default_imagescXlim)

parse(p, filePath, fileName, varName, vargs{:})


end
