function p = parse_input_getGlobalUnitIdx( filePath, fileName, eventName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function

%parse input, and extract name-value pairs for the main function
% 'runPCA'
default_frHighPass  = 1;  % low FR cut (3Hz)
default_frLowPass   = 60; % high FR cut (60Hz) % low pass might not be necessary
default_xcorThresholdPer = 0.2; % default cross-correlation spike cooccurrence tolerance (i.e., 20% of the total spike count of the fewer spiking unit of the pair)
default_redefineXcorrUnits = false;

p = inputParser; % create parser object
addRequired(p,'filePath'); % file directory
addRequired(p,'fileName'); % file name
addRequired(p,'eventName'); % event name, e.g. reach, tagLaser - the events around which the psths were taken

addParameter(p,'frHighPass', default_frHighPass)
addParameter(p,'frLowPass', default_frLowPass) % FR high cut might not be necessary
addParameter(p,'xcorThresholdPer', default_xcorThresholdPer)
addParameter(p,'redefineXcorrUnits', default_redefineXcorrUnits)

parse(p, filePath, fileName, eventName, vargs{:})
end