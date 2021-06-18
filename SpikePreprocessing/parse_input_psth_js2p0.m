function p = parse_input_psth_js2p0( filePath, fileInfo, probeDepth, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
%parse input, and extract name-value pairs for the main function
% 'PSTH_rasters'

default_probeAngle = 10;       % default probe angle in degree
default_numbSiteProbe = 384;   % default number of sites per probe(s) (can be an arrary, if there are multiple probes used
default_psthPlotFlag = false;  % default logic indicating psth draw or not
default_reachWin = [3e3 2e3];  % default time window for reach psth
default_rewardWin = [3e3 2e3]; % default time window for reward psth
default_tagLaserWin = [5e3 5e3]; % default time window for tagLaser psth
default_probeType = 'imec';    % default probe type imec
default_laserUsed = 'true';    % default laser used logic

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'fileInfo');
addRequired(p,'probeDepth');
addParameter(p,'probeAngle', default_probeAngle)
addParameter(p,'numbSiteProbe', default_numbSiteProbe)
addParameter(p,'psthPlotFlag', default_psthPlotFlag)
addParameter(p,'reachWin', default_reachWin)
addParameter(p,'rewardWin', default_rewardWin)
addParameter(p,'tagLaserWin', default_tagLaserWin)
addParameter(p,'probeType', default_probeType)
addParameter(p,'laserUsed', default_laserUsed)

parse(p,filePath, fileInfo, probeDepth, vargs{:})
end