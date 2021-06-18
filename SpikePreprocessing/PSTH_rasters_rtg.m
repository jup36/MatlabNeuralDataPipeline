function PSTH_rasters_rtg( filePath, fileInfo, probeDepth, varargin )
%PSTH_rasters takes behavioral timestamps andvgenerates psths aligned to 
% each event and saves the outcome psths in the filePath. 

p = parse_input_psth_rtg(filePath, fileInfo, probeDepth, varargin); % parse input
%filePath = '/Volumes/Beefcake/Junchol_Data/jayReachToGrasp/M314Sim1GtACr2Junchol/M314_20200427_1000um_g0'; 
%fileInfo = 'M314_20200427_1000um'; 
%p = parse_input_psth_rtg(filePath,fileInfo,1000,{'probeAngle',0}) % when running line-by-line

%% Load files 
cd(p.Results.filePath)   % change directory to the data folder

% get behavioral data
behFile = dir(fullfile(p.Results.filePath,'BehVariablesJayRtg.mat')); 
load(fullfile(behFile.folder,behFile.name),'ts');

% get meta 
disp('Select the meta file!!')
[metaFileSel,metaPathSel] = uigetfile('*.meta',p.Results.filePath); 
meta = ReadMeta(metaFileSel, metaPathSel); % read out the meta file
% get geometry using meta
%SGLXMetaToCoords(meta, metaFileSel) % make a chanMap file from meta, set outType=1 for ks2 format
if ispc
    load(fullfile('S:\Junchol_Data\jayReachToGrasp\jay_64_16_4.mat'),'xcoords','ycoords'); 
elseif ismac
    load(fullfile('/Volumes/Beefcake/Junchol_Data/jayReachToGrasp/jay_64_16_4.mat'),'xcoords','ycoords'); 
end

%load(fullfile('/Volumes/Beefcake/Junchol_Data/jayReachToGrasp/jay_64_16_4.mat')); 
geometry = [xcoords, ycoords]; % probe x, y coordinates 

%% kilosort-phy
spike_times = double(readNPY(fullfile(filePath,'spike_times.npy'))); % timestamp of all spikes, (n_spike, 1)
spike_template = double(readNPY(fullfile(filePath,'spike_templates.npy'))); % automated cluster of all spikes, (n_spike, 1)
spike_clusters = double(readNPY(fullfile(filePath, 'spike_clusters.npy')));  % final cluster of all spikes, (n_spike, 1)
template = readNPY(fullfile(filePath, 'templates.npy')); % template waveform of all clusters, (n_original_cluster, n_timepoint, all_valid_channel, i.e., 64)
channel_map = readNPY(fullfile(filePath, 'channel_map.npy')); % maps valid sites to actual sites

% find the main channel of each template which has the greatest abs amplitude
[~, mainSiteTemplate] = max(max(abs(template), [], 2), [], 3);
actualSiteTemplate = channel_map(mainSiteTemplate)+1; 

% load cluster data (final cluster IDs after the manual curation)
dc = dir('cluster_group*'); 
if contains(dc.name,'.tsv')
     movefile('cluster_group.tsv','cluster_group.csv') 
end

try
    cnTab = readtable(fullfile(dc.folder,'cluster_group.csv')); 
    cn = table2cell(cnTab); 
    inPhy = cell2mat(cellfun(@(a) strcmpi(a,'good'), cn(:,2), 'un',0));
    unitNumber = cell2mat(cn(inPhy,1)); % final good unit number list
    nP = length(unitNumber);
catch
    cn_name = fullfile(filePath, 'cluster_groups.csv');
    fid = fopen(cn_name, 'r');
    cn = textscan(fid, '%f%s%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'Headerlines', 1, 'EndOfLine', '\r\n');
    fclose(fid);
    % pick only good units
    inPhy = strcmp(cn{2}, 'good'); 
    unitNumber = cn{1}(inPhy); % final good unit number list
    nP = length(unitNumber);    
end

dvCosConvert   = cos(p.Results.probeAngle/180*pi);  % if probe was angled, probe coordinates need to be corrected 

if length(p.Results.probeDepth)~=length(p.Results.numbSiteProbe)
    error('Check the variables; probeDepth & numbChProbe!')
end

% Assign probe id to channels (as multiple probes might be used)
whichProbe = zeros(sum(p.Results.numbSiteProbe),1); % probe identity for each channel  

probeId = 1; % one-based (first probe gets 1)
for s = 1:sum(p.Results.numbSiteProbe) % increment all probe channels
    if s <= sum(p.Results.numbSiteProbe(1:probeId))
        whichProbe(s,1) = probeId;
    elseif s > sum(p.Results.numbSiteProbe(1:probeId))
        whichProbe(s,1) = probeId+1;
        probeId = probeId+1;
    end
end
clearvars s 

%% Process spike times data and generate PSTH and Rasters
maxST = max(spike_times)/(str2double(meta.niSampRate)/1000); % in ms
spkTimes = struct; % the structure to contain spike times
for u = 1:length(unitNumber) % increment valid clusters (units)
    spkIdx = unitNumber(u);  %S_clu.cviSpk_clu{u}; % spikeIDs of the current cluster 
    templateIdx = mode(spike_template(spike_clusters==spkIdx))+1; % template Id
    if strcmp(meta.typeThis, 'nidq')
        spkTimes(u).spkTimes = spike_times(spike_clusters==spkIdx)/(str2double(meta.niSampRate)/1000); % divided by the sampling rate: 25kHz (str2num(meta.niSampRate)/1000)
    end
    
    spkTimes(u).origClusId  = spkIdx; % original cluster ID back in kilosort/phy
    spkTimes(u).maxSite     = actualSiteTemplate(templateIdx);  % assign the current cluster to a site of the probe
    spkTimes(u).geometry(1) = geometry(spkTimes(u).maxSite,1); % horizontal geometry
    spkTimes(u).geometry(2) = p.Results.probeDepth(whichProbe(spkTimes(u).maxSite))-geometry(spkTimes(u).maxSite,2); % vertical geometry
    spkTimes(u).geometry(2) = spkTimes(u).geometry(2).*dvCosConvert; % corrected dv (original dv multiplied by cosine(probe angle))
    spkTimes(u).isStr = spkTimes(u).geometry(2)>2000;  % logical for striatum 
    spkTimes(u).template = template(templateIdx,:,mainSiteTemplate(templateIdx)); % the spike template for each cluster
    spkTimes(u).templateId = templateIdx; 
end
clearvars u 

if isfield(spkTimes,'spkTimesCr')
  spkTimes = rmfield(spkTimes,'spkTimes'); % drop the uncorrected time points
end
spkTimesCell = struct2cell(spkTimes'); % the entire spike times converted into a cell 

spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(3,:))<=64); % the CTX spike times cell (1st probe)
%spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(3,:))>64);  % the STR spike times cell (2nd probe)

%spkTimesCellStrCtx = [spkTimesCellSTR,  spkTimesCellCTX];

%% get psths
% binned spike count CTX
<<<<<<< HEAD
cueRangeC = arrayfun(@(a) a-2000:a+2000, ts.cue, 'un', 0); 
cueNoLaserI = cell2mat(cellfun(@(a) sum(ismember(ts.laserCue2s,a))==0, cueRangeC, 'un', 0)); 
ts.cueNoLaser = ts.cue(cueNoLaserI); 

if ~isempty(ts.tagLaser1s)
    tagLaser1s     = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, ts.tagLaser1s', ts.cue'-1000, 1, [5e3 5e3], -1, p.Results.psthPlotFlag );
=======
cueRangeC = arrayfun(@(a) a-3000:a+3000, ts.cue, 'un', 0); 
cueNoLaserI = cell2mat(cellfun(@(a) sum(ismember([ts.laserCue2s, ts.laserOnly2s, ts.tagLaser1s],a))==0, cueRangeC, 'un', 0)); 
ts.cueNoLaser = ts.cue(cueNoLaserI); 

if ~isempty(ts.tagLaser1s)
    tagLaser1s     = psthBINcellRmvArtifactNearZero( p.Results.fileInfo, 'M1', spkTimesCellCTX, ts.tagLaser1s', ts.cue'-1000, 1, [5e3 5e3], -1, p.Results.psthPlotFlag );
>>>>>>> master
    binSpkCountCTX.tagLaser1s = tagLaser1s;
end

if  ~isempty(ts.laserCue2s)
    laserCue2s     = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, ts.laserCue2s', ts.cue'-1000, 1, [5e3 5e3], -1, p.Results.psthPlotFlag );
    binSpkCountCTX.laserCue2s  = laserCue2s;
end

if  ~isempty(ts.laserOnly2s)
<<<<<<< HEAD
    laserOnly2s    = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, ts.laserOnly2s', ts.cue'-1000, 1, [5e3 5e3], -1, p.Results.psthPlotFlag );
=======
    laserOnly2s    = psthBINcellRmvArtifactNearZero( p.Results.fileInfo, 'M1', spkTimesCellCTX, ts.laserOnly2s', ts.cue'-1000, 1, [5e3 5e3], -1, p.Results.psthPlotFlag );
>>>>>>> master
    binSpkCountCTX.laserOnly2s = laserOnly2s;
end

if  ~isempty(ts.cueNoLaser)
    cueNoLaser    = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, ts.cueNoLaser', ts.cue'-1000, 1, [5e3 5e3], -1, p.Results.psthPlotFlag );
    binSpkCountCTX.cueNoLaser = cueNoLaser;
end

saveName = strcat('binSpkCountCTX',p.Results.fileInfo);
save(fullfile(p.Results.filePath,saveName),'-struct','binSpkCountCTX') % save the fields of the structure separately 
save(fullfile(p.Results.filePath,saveName), 'ts', '-append') % append the behavioral timestamps

%% Individual unit raster plot 
<<<<<<< HEAD
unit = 6; % 32, 37, 40, 45, 66, 69 (M314_20200427_1000um)
spikeRasterGramm( [5e3 5e3], {'tagLaser1s'}, [3e3 3e3], [tagLaser1s.SpkTimes{unit};laserOnly2s.SpkTimes{unit}]);
spikeRasterGramm( [5e3 5e3], {'tagLaser1s'}, [3e3 3e3], tagLaser1s.SpkTimes{unit});
spikeRasterGramm( [5e3 5e3], {'laserOnly2s'}, [3e3 3e3], laserOnly2s.SpkTimes{unit});
%spikeRasterGramm( [5e3 5e3], {'tagLaser1s'}, [3e3 3e3], binSpkCountCTX.tagLaser1s.SpkTimes{unit});
%spikeRasterGramm( [5e3 5e3], {'laserOnly2s'}, [3e3 3e3], binSpkCountCTX.laserOnly2s.SpkTimes{unit});
%print( fullfile(filePath,'Figure',strcat(fileInfo,'_',sprintf('unit#%d',unit),'tagLaser1sLaserOnly2s')), '-dpdf','-painters', '-bestfit')
=======
%unit = 1; % 32, 37, 40, 45, 66, 69 (M314_20200427_1000um)
%spikeRasterGramm( [5e3 5e3], {'tagLaser1s'}, [3e3 3e3], [tagLaser1s.SpkTimes{unit};laserOnly2s.SpkTimes{unit}]);
%spikeRasterGramm( [5e3 5e3], {'tagLaser1s'}, [3e3 3e3], tagLaser1s.SpkTimes{unit});
%spikeRasterGramm( [5e3 5e3], {'laserOnly2s'}, [3e3 3e3], laserOnly2s.SpkTimes{unit});
%print( fullfile(filePath,'Figure',strcat(fileInfo,'_',sprintf('unit#%d',unit),'LaserOnly2s')), '-dpdf','-painters', '-bestfit')
>>>>>>> master
%unit = unit+1;

%unit = 139; % str unit 128 (laser activated)
%spikeRasterGramm( [1e3 3e3], {'reach','stim','stimReach'}, binSpkCountSTRReach(unit).SpkTimes, binSpkCountSTRstmLaser(unit).SpkTimes,binSpkCountSTRstmReach(unit).SpkTimes );

<<<<<<< HEAD
%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
>>>>>>> master
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
end


%%%%%%%%%%%%%%%%%%%%%%%
%  NON-NESTED FUCNTIONS
%%%%%%%%%%%%%%%%%%%%%%%
function [ geometry ] = getHH3geom
%This function returns the geometry of the imec3 option3 probe. 

% Order of the probe sites in the recording file
geometry = zeros(64,2);
chanMap = [49,50,51,52,54,53,56,55,58,57,60,...
           59,62,61,64,63,16,15,14,13,11,12,...
           9,7,5,3,1,10,6,4,2,8,25,31,29,27,...
           23,32,30,28,26,24,21,22,20,19,18,...
           17,34,33,36,35,38,37,40,39,42,41,...
           44,43,45,46,47,48]; % 1-indexed
ycoords = 0:20:20*63; 
       
% Site location in micrometers (x and y)
geometry(chanMap,2) = ycoords;
end

function [ meta ] = getmeta
%This helper function returns meta data 
% using ReadMeta (by Bill Karsh) that can read out the meta file generated by spikeGLX
% ReadMeta requires 1) name of the bin file, 2) file directory of the
% binfile (given as the current file directory - pwd)

fileList = dir('*ap.bin'); % identify the *.ap.bin file 

if isempty(fileList)
    fileList = dir('*nidq.bin');
end

if length(fileList)==1
    meta = ReadMeta(fileList.name, pwd); % read out the meta file 
elseif length(fileList)>1 
    error('Make sure there is only one *ap.bin file in the current folder')
elseif isempty(fileList)
    error('No *ap.bin file was detected!')
end

end

function [ S_clu,viTime_spk,viSite_spk ] = getjrcmatVar
%This function loads the 'S_Clu' (the structure containing spike
% information), and 'viTime_spk' (timestamp per spike - to be divided by the sampling rate for conversion to actual time).  
% 'viSite_spk' (site with the peak spike amplitude).
% in the '*jrc.mat' file.

fileList = dir('*_jrc.mat'); % identify the *.ap.bin file 

if length(fileList)==1 
    load(fileList.name,'S_clu','viTime_spk','viSite_spk');
else 
    error('Check the *jrc.mat file; one *_jrc.mat file must exist in the working folder!')
end

end