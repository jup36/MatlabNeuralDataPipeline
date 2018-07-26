function p = parse_input_stimEffect( filePath, fileName, probeDepth, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
%parse input, and extract name-value pairs for the main function
% 'stimEffect.m'

default_reachWinEdges = -2e3+1:2e3; % default reach window
default_reachOnTime   = 0;   % time at which reach is on
default_reachDur      = 750; % default reach duration

default_tagWinEdges = -5e3+1:5e3; % default tag window
default_tagOnTime   = 0;   % time at which tag stim is on
default_tagDur      = 1000; % default tag stim duration
default_tagColorAxis = [-5 5]; % default tag color axis
default_binStepSize = 10;  % 10 ms to be used for spike counts
default_binSize = 50; % 50ms bin to be used for spike counts

default_lowFRcut  = 0.5;    % to exclude units with low FR
default_colorAxis = [-5 5]; % color axis for stim effect probe colormap

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'fileName');
addRequired(p,'probeDepth');

addParameter(p,'reachWinEdges', default_reachWinEdges)
addParameter(p,'reachOnTime', default_reachOnTime)
addParameter(p,'reachDur', default_reachDur)
addParameter(p,'tagWinEdges', default_tagWinEdges)
addParameter(p,'tagOnTime', default_tagOnTime)
addParameter(p,'tagDur', default_tagDur)
addParameter(p,'lowFRcut', default_lowFRcut)
addParameter(p,'colorAxis', default_colorAxis)
addParameter(p,'tagColorAxis', default_tagColorAxis)
addParameter(p,'binStepSize', default_binStepSize)
addParameter(p,'binSize', default_binSize)

parse(p, filePath, fileName, probeDepth, vargs{:})

end