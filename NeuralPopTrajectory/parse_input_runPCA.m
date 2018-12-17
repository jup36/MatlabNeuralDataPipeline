function p = parse_input_runPCA( filePath, fileName, saveNameTag, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function

%parse input, and extract name-value pairs for the main function
% 'runPCA'

default_frHighPass  = 1;  % low FR cut (3Hz)
default_frLowPass   = 60; % high FR cut (60Hz) % low pass might not be necessary
default_timeRange   = [-0.5e3 0.5e3]; % in ms (e.g. -500 to 500 ms relative to the event at time 0)
default_pcaBinSize = 50;  % bin size to be used for PCA (25ms)
default_pcaDim     = 5;   % default dimensionality for PCA
default_pcaKernSD  = 50;  % default SD for kernel smoothing in PCA
default_fullCrossVal= false; % default boolean to run full cross validation or not
default_xcorThresholdPer = 0.2; % default cross-correlation spike cooccurrence tolerance (i.e., 20% of the total spike count of the fewer spiking unit of the pair)
default_crossValFolds = 4; % the number of cross-validation fold
default_baseSubtrt = true; % the logical to subtract the baseline activity (spike count) or not (by default, subtract the mean spike counts across the bins corresponding to the leftmost 1-s window)
default_sqrtSC = false;    % the logical to take sqrt of the spike counts or not before conducting dimensionality reduction
default_redefineXcorrUnits = false;

p = inputParser; % create parser object
addRequired(p,'filePath'); % file directory
addRequired(p,'fileName'); % file name to open
addRequired(p,'saveNameTag'); % save name to save pcaResult

addParameter(p,'frHighPass', default_frHighPass)
addParameter(p,'frLowPass', default_frLowPass) % FR high cut might not be necessary
addParameter(p,'timeRange', default_timeRange)
addParameter(p,'pcaBinSize', default_pcaBinSize)
addParameter(p,'pcaDim', default_pcaDim)
addParameter(p,'pcaKernSD', default_pcaKernSD)
addParameter(p,'fullCrossVal', default_fullCrossVal)
addParameter(p,'xcorThresholdPer', default_xcorThresholdPer)
addParameter(p,'crossValFolds', default_crossValFolds)
addParameter(p,'baseSubtrt', default_baseSubtrt)
addParameter(p,'sqrtSC', default_sqrtSC)
addParameter(p,'redefineXcorrUnits', default_redefineXcorrUnits)

parse(p, filePath, fileName, saveNameTag, vargs{:})

end