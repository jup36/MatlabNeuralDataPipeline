function p = parse_input_stimEffectReach( filePath, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
% parse input, and extract name-value pairs

default_useAllTrials = true;  % logical to indicate using all trials in analyses
default_selectTrials = 1:150; % default select trials
default_minReachDurCut = 50;  % default minimum reach duration cut off (e.g. 50 ms) - to deselect unless the detected reach duration is greater than this value
default_reachBeforeLastReward = true; % detect only the reaches occur before the last reward delivery

p = inputParser; % create parser object

addRequired(p,'filePath');
addParameter(p,'useAllTrials', default_useAllTrials)
addParameter(p,'selectTrials', default_selectTrials)
addParameter(p,'minReachDurCut', default_minReachDurCut)
addParameter(p,'reachBeforeLastReward', default_reachBeforeLastReward)

parse(p,filePath,vargs{:})

end
