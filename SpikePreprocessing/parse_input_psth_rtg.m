function p = parse_input_psth_rtg( filePath, fileInfo, probeDepth, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
%parse input, and extract name-value pairs for the main function
% 'PSTH_rasters'

default_probeAngle = 0;       % default probe angle in degree
default_numbSiteProbe = 64;   % default number of sites per probe(s) (can be an arrary, if there are multiple probes used
default_psthPlotFlag = false; % default logic indicating psth draw or not
default_cueWin = [3e3 3e3];   % default time window for cue onset
default_tagLaser1sWin = [5e3 5e3]; % default time window for tagLaser (1s) psth
default_laserOnly2sWin = [5e3 5e3]; % default time window for laserOnly (2s) psth
default_laserCue2sWin = [5e3 5e3];  % default time window for laserCue (2s) psth

default_probeType = '64_16_4';  % default probe type

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'fileInfo');
addRequired(p,'probeDepth');
addParameter(p,'probeAngle', default_probeAngle)
addParameter(p,'numbSiteProbe', default_numbSiteProbe)
addParameter(p,'psthPlotFlag', default_psthPlotFlag)
addParameter(p,'cueWin', default_cueWin)
addParameter(p,'tagLaser1sWin', default_tagLaser1sWin)
addParameter(p,'laserOnly2sWin', default_laserOnly2sWin)
addParameter(p,'laserCue2sWin', default_laserCue2sWin)
addParameter(p,'default_probeType', default_probeType)

parse(p,filePath, fileInfo, probeDepth, vargs{:})
end