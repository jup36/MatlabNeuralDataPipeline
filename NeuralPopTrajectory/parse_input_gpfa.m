function p = parse_input_gpfa( filePath, fileName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function

%parse input, and extract name-value pairs for the main function
% 'PSTH_rasters'

default_frHighPass  = 2;  % low FR cut (2Hz)
default_frLowPass   = 50; % high FR cut (50Hz)
default_timeRange   = [-0.5e3 0.5e3]; % in ms (e.g. -500 to 500 ms relative to the event at time 0)
default_gpfaBinSize = 50; % bin size to be used for GPFA
default_gpfaDim     = 5;  % default dimensionality for GPFA
default_kernSD      = 50; % default SD for kernel smoothing in GPFA
default_fullCrossVal= false; % default boolean to run full cross validation or not

p = inputParser; % create parser object
addRequired(p,'filePath'); % file directory
addRequired(p,'fileName'); % file name to open

addParameter(p,'frHighPass', default_frHighPass)
addParameter(p,'frLowPass', default_frLowPass)
addParameter(p,'timeRange', default_timeRange)
addParameter(p,'gpfaBinSize', default_gpfaBinSize)
addParameter(p,'gpfaDim', default_gpfaDim)
addParameter(p,'kernSD', default_kernSD)
addParameter(p,'fullCrossVal', default_fullCrossVal)

parse(p, filePath, fileName, vargs{:})

end


