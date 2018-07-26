function p = parse_input_runGPFA( filePath, fileName, saveNameTag, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function

%parse input, and extract name-value pairs for the main function
% 'runGPFA'

default_frHighPass  = 3;  % low FR cut (3Hz)
default_frLowPass   = 60; % high FR cut (60Hz) % low pass might not be necessary
default_timeRange   = [-0.5e3 0.5e3]; % in ms (e.g. -500 to 500 ms relative to the event at time 0)
default_gpfaBinSize = 50; % bin size to be used for GPFA (25ms)
default_gpfaDim     = 5;  % default dimensionality for GPFA
default_gpfaKernSD  = 50; % default SD for kernel smoothing in GPFA
default_fullCrossVal= false; % default boolean to run full cross validation or not
default_xcorThresholdPer = 0.2; % default cross-correlation spike cooccurrence tolerance (i.e., 20% of the total spike count of the fewer spiking unit of the pair)
default_crossValFolds = 4;  % the number of cross-validation fold

p = inputParser; % create parser object
addRequired(p,'filePath'); % file directory
addRequired(p,'fileName'); % file name to open
addRequired(p,'saveNameTag'); % save name to save gpfaResult

addParameter(p,'frHighPass', default_frHighPass)
addParameter(p,'frLowPass', default_frLowPass) % FR high cut might not be necessary
addParameter(p,'timeRange', default_timeRange)
addParameter(p,'gpfaBinSize', default_gpfaBinSize)
addParameter(p,'gpfaDim', default_gpfaDim)
addParameter(p,'gpfaKernSD', default_gpfaKernSD)
addParameter(p,'fullCrossVal', default_fullCrossVal)
addParameter(p,'xcorThresholdPer', default_xcorThresholdPer)
addParameter(p,'crossValFolds', default_crossValFolds)

parse(p, filePath, fileName, saveNameTag, vargs{:})

end


