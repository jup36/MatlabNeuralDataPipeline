function PSTH_rasters_js2p0( filePath, fileInfo, probeDepth, varargin )
%PSTH_rasters takes behavioral timestamps stored in BehVariables.mat and
% generates psths aligned to each event and saves the outcome psths in the
% filePath. To run PSTH_rasters 'behaviorTimestamps.m' must be run first
% and its outcome 'BehVariables.mat' must exist in the filePath. 
% Modified on 6/18/18 to save the combined (all imec channels)
% 'binSpkCountStrCtx*.mat'
% Modified in Jan/19 to change the baseline period for tagLaser PSTHs from
% reachStart to tagLaser. 

p = parse_input_psth_js2p0(filePath, fileInfo, probeDepth, varargin); % parse input
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/imec'; 
%fileInfo = 'WR40_081919'; 
%p = parse_input_psth_js2p0(filePath,'WR40_081919',4200,{'probeAngle',10,'numbSiteProbe',384,'psthPlotFlag',false,'reachWin',[3e3 3e3],'rewardWin',[3e3 3e3],'tagLaserWin',[5e3 5e3]}) % when running line-by-line

%% Load files 
cd(p.Results.filePath)   % change directory to the data folder

% get behavioral data
behFile = dir(fullfile(p.Results.filePath,'BehVariablesJs.mat')); 
if length(behFile)>1 || isempty(behFile)
    disp('Select BehVariablesJs.mat!')
    [behFileSelect,behPathSelect] = uigetfile('*.mat',p.Results.filePath);
    evt = load(fullfile(behPathSelect,behFileSelect),'evtIdx1k');
    evt = evt.('evtIdx1k'); 
    clearvars evtIdx1k
else
    evt = load(fullfile(behFile.folder,behFile.name),'evtIdx1k');
    evt = evt.('evtIdx1k'); 
    clearvars evtIdx1k
end

jkvtFile = dir(fullfile(p.Results.filePath,'jsTime1k_KinematicsTrajectories.mat')); 
if length(jkvtFile)>1 || isempty(jkvtFile)
    disp('Select jsTime1k_KinematicsTrajectories.mat!')
    [jkvtFileSelect,jkvtPathSelect] = uigetfile('.mat',p.Results.filePath);
    jkvt = load(fullfile(jkvtPathSelect,jkvtFileSelect),'jkvt');
    jkvt = jkvt.('jkvt'); 
else
    jkvt = load(fullfile(jkvtFile.folder,jkvtFile.name),'jkvt');
    jkvt = jkvt.('jkvt'); 
end

% get meta 
disp('Select the meta file!!')
[metaFileSel,metaPathSel] = uigetfile('*.meta',p.Results.filePath); 
meta = ReadMeta(metaFileSel, metaPathSel); % read out the meta file
% get geometry using meta
SGLXMetaToCoords(meta, metaFileSel) % make a chanMap file from meta, set outType=1 for ks2 format
geomFile = dir(fullfile(p.Results.filePath,'*_kilosortChanMap.mat'));  % look for '*_kilosortChanMap.mat' file 
load(fullfile(geomFile.folder,geomFile.name),'xcoords','ycoords'); 
geometry = [xcoords, ycoords]; % probe x, y coordinates 

%% kilosort-phy
phyPath = dir(fullfile(p.Results.filePath,'*/*spike_times.npy')); 
if isempty(phyPath)
    disp('Select a folder containing spike sorting results with Phy!!')
    phyPath = uigetdir(p.Results.filePath); 
elseif size(phyPath,1)>1
    phyPath = phyPath(cellfun(@(a) contains(a,p.Results.probeType(1:2),'IgnoreCase',true), {phyPath(:).folder})).folder; 
end

spike_times = double(readNPY(fullfile(phyPath, 'spike_times.npy'))); % timestamp of all spikes, (n_spike, 1)
spike_template = double(readNPY(fullfile(phyPath, 'spike_templates.npy'))); % automated cluster of all spikes, (n_spike, 1)
spike_clusters = double(readNPY(fullfile(phyPath, 'spike_clusters.npy')));  % final cluster of all spikes, (n_spike, 1)
template = readNPY(fullfile(phyPath, 'templates.npy')); % template waveform of all clusters, (n_original_cluster, n_timepoint, all_valid_channel, i.e., 64)
channel_map = readNPY(fullfile(phyPath, 'channel_map.npy')); % maps valid sites to actual sites

% find the main channel of each template which has the greatest abs amplitude
[~, mainSiteTemplate] = max(max(abs(template), [], 2), [], 3);
actualSiteTemplate = channel_map(mainSiteTemplate)+1; 

% load cluster data (final cluster IDs after the manual curation)
cn_name = fullfile(phyPath, 'cluster_groups.csv');
fid = fopen(cn_name, 'r');
cn = textscan(fid, '%f%s%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'Headerlines', 1, 'EndOfLine', '\r\n');
fclose(fid);

% pick only good units
inPhy = strcmp(cn{2}, 'good'); 
unitNumber = cn{1}(inPhy); % final good unit number list
nP = length(unitNumber);

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
sTcr = [jkvt(:).trEndImec]-[jkvt(:).trEnd]; % to correct the spike times due to asynchrony between the two recording systems
maxST = max(spike_times)/(str2double(meta.imSampRate)/1000); 
spkTimes = struct; % the structure to contain spike times
for u = 1:length(unitNumber) % increment valid clusters (units)
    spkIdx = unitNumber(u);  %S_clu.cviSpk_clu{u}; % spikeIDs of the current cluster 
    templateIdx = mode(spike_template(spike_clusters==spkIdx))+1; % template Id
    if strcmp(meta.typeThis, 'imec')
        spkTimes(u).spkTimes = spike_times(spike_clusters==spkIdx)/(str2double(meta.imSampRate)/1000); % divided by the sampling rate: 30kHz (str2num(meta.imSampRate)/1000)
        [~,~,trBin] = histcounts(spkTimes(u).spkTimes, [1 jkvt(:).trEnd maxST+1]); 
        spkTimes(u).spkTimesCr = spkTimes(u).spkTimes-sTcr(min(trBin,length(sTcr)))'; 
        spkTimes(u).spkTimesCr = max(1, spkTimes(u).spkTimesCr); 
    elseif strcmp(meta.typeThis, 'nidq')
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
if strcmp(meta.typeThis, 'imec')          % in case recording via imec  
    spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(5,:))==0); % the CTX spike times cell (cells from sites above the strCtxBorder)
    spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(5,:))==1); % the STR spike times cell (cells from sites below the strCtxBorder )
elseif strcmp(meta.typeThis, 'nidq')      % in case recording from one or two apig probes
    spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(3,:))<=64); % the CTX spike times cell (1st probe)
    spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(3,:))>64);  % the STR spike times cell (2nd probe)
end

spkTimesCellStrCtx = [spkTimesCellSTR,  spkTimesCellCTX]; 

%% Specify all the time points to align neural data onto 
% define pull Starts/Stops
evt.trJsReady = [jkvt.trJsReady]'; % joystick ready timePoint per trial
pullTrsIdx = cellfun(@(c)strcmpi(c,'sp'), {jkvt.trialType}); % successfull pull trials

for t = 1:length(jkvt)
    % detect pull start/stop of a successful trials
    if pullTrsIdx(t) && ~isempty(jkvt(t).movKins.pullStart) && ~isempty(jkvt(t).movKins.pullStop)
        jkvt(t).pullStarts = jkvt(t).trJsReady + jkvt(t).movKins.pullStart; 
        jkvt(t).pullStops  = jkvt(t).trJsReady + jkvt(t).movKins.pullStop; 
        % detect reach start (hand lift) of each successful trial 
        tmphTrjRstartT = jkvt(t).vFrameTime(jkvt(t).hTrjRstart); 
        jkvt(t).rStartToPull = tmphTrjRstartT(find(tmphTrjRstartT<jkvt(t).pullStarts,1,'last')); 
        % detect reach stop 
        tmphTrjRstopT = jkvt(t).vFrameTime(jkvt(t).hTrjRstop); 
        jkvt(t).rStopToPull = tmphTrjRstopT(find(tmphTrjRstopT>jkvt(t).pullStarts,1,'last')); 
    end
end
clearvars t

% binned spike count CTX 
pullStarts     = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, [jkvt.pullStarts]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag ); 
pullStarts.trI = find(~cellfun(@isempty, {jkvt.pullStarts}));  
pullStops      = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, [jkvt.pullStops]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
pullStops.trI  = find(~cellfun(@isempty, {jkvt.pullStops}));  
rStartToPull   = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, [jkvt.rStartToPull]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
rStartToPull.trI = find(~cellfun(@isempty, {jkvt.rStartToPull}));  
rStopToPull    = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, [jkvt.rStopToPull]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
rStopToPull.trI = find(~cellfun(@isempty, {jkvt.rStopToPull}));  
reward         = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, evt.rwdIdx', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );     % entire rewardDelivery
reward.trI     = find([jkvt.rewarded]);  
stmLaser       = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, evt.stimLaserRiseIdx', evt.stimLaserRiseIdx', 1, [1e3 5e3], -1, p.Results.psthPlotFlag ); % laser stim trials
tagLaser       = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, evt.tagLaserRiseIdx', evt.tagLaserRiseIdx', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag ); % laser tag trials

binSpkCountCTX.pullStarts = pullStarts; 
binSpkCountCTX.pullStops  = pullStops; 
binSpkCountCTX.rStartToPull = rStartToPull; 
binSpkCountCTX.rStopToPull = rStopToPull; 
binSpkCountCTX.reward = reward; 
binSpkCountCTX.stmLaser = stmLaser; 
binSpkCountCTX.tagLaser = tagLaser; 
binSpkCountCTX.meta = meta; 
binSpkCountCTX.p = p; 
binSpkCountCTX.spkTimesCell = spkTimesCell; % just to save the cell

saveName = strcat('binSpkCountCTX',p.Results.fileInfo);
save(fullfile(p.Results.filePath,saveName),'-struct','binSpkCountCTX') % save the fields of the structure separately 
save(fullfile(p.Results.filePath,saveName), 'evt', 'jkvt', '-append') % append the behavioral timestamps
clearvars saveName binSpkCountCTX pullStarts pullStops trStart trEnd reward stmLaser tagLaser rStartToPull rStopToPull 

%% binned spike count STR 
pullStarts     = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, [jkvt.pullStarts]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag ); 
pullStarts.trI = find(~cellfun(@isempty, {jkvt.pullStarts}));  
pullStops      = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, [jkvt.pullStops]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
pullStops.trI  = find(~cellfun(@isempty, {jkvt.pullStops}));  
rStartToPull   = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, [jkvt.rStartToPull]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
rStartToPull.trI = find(~cellfun(@isempty, {jkvt.rStartToPull}));  
rStopToPull    = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, [jkvt.rStopToPull]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
rStopToPull.trI = find(~cellfun(@isempty, {jkvt.rStopToPull}));  
reward         = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, evt.rwdIdx', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );     % entire rewardDelivery
reward.trI     = find([jkvt.rewarded]);  
stmLaser       = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, evt.stimLaserRiseIdx', evt.stimLaserRiseIdx', 1, [1e3 5e3], -1, p.Results.psthPlotFlag ); % laser stim trials
tagLaser       = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, evt.tagLaserRiseIdx', evt.tagLaserRiseIdx', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag ); % laser tag trials

binSpkCountSTR.pullStarts = pullStarts; 
binSpkCountSTR.pullStops  = pullStops; 
binSpkCountSTR.rStartToPull = rStartToPull; 
binSpkCountSTR.rStopToPull = rStopToPull; 
binSpkCountSTR.reward = reward; 
binSpkCountSTR.stmLaser = stmLaser; 
binSpkCountSTR.tagLaser = tagLaser; 
binSpkCountSTR.meta = meta; 
binSpkCountSTR.p = p; 
binSpkCountSTR.spkTimesCell = spkTimesCell; % just to save the cell

saveName = strcat('binSpkCountSTR',p.Results.fileInfo);
save(fullfile(p.Results.filePath,saveName),'-struct','binSpkCountSTR') % save the fields of the structure separately 
save(fullfile(p.Results.filePath,saveName), 'evt', 'jkvt', '-append') % append the behavioral timestamps
clearvars saveName binSpkCountSTR pullStarts pullStops trStart trEnd reward stmLaser tagLaser rStartToPull rStopToPull 

% binned spike count STRCTX 
pullStarts     = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, [jkvt.pullStarts]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag ); 
pullStarts.trI = find(~cellfun(@isempty, {jkvt.pullStarts}));  
pullStops      = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, [jkvt.pullStops]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
pullStops.trI  = find(~cellfun(@isempty, {jkvt.pullStops}));  
rStartToPull   = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, [jkvt.rStartToPull]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
rStartToPull.trI = find(~cellfun(@isempty, {jkvt.rStartToPull}));  
rStopToPull    = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, [jkvt.rStopToPull]', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );
rStopToPull.trI = find(~cellfun(@isempty, {jkvt.rStopToPull}));  
reward         = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, evt.rwdIdx', evt.trJsReady-1000, 1, [3e3 2e3], -1, p.Results.psthPlotFlag );     % entire rewardDelivery
reward.trI     = find([jkvt.rewarded]);  
stmLaser       = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, evt.stimLaserRiseIdx', evt.stimLaserRiseIdx', 1, [1e3 5e3], -1, p.Results.psthPlotFlag ); % laser stim trials
tagLaser       = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, evt.tagLaserRiseIdx', evt.tagLaserRiseIdx', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag ); % laser tag trials

binSpkCountSTRCTX.pullStarts = pullStarts; 
binSpkCountSTRCTX.pullStops  = pullStops; 
binSpkCountSTRCTX.rStartToPull = rStartToPull; 
binSpkCountSTRCTX.rStopToPull = rStopToPull; 
binSpkCountSTRCTX.reward = reward; 
binSpkCountSTRCTX.stmLaser = stmLaser; 
binSpkCountSTRCTX.tagLaser = tagLaser; 
binSpkCountSTRCTX.meta = meta; 
binSpkCountSTRCTX.p = p; 
binSpkCountSTRCTX.spkTimesCell = spkTimesCell; % just to save the cell

saveName = strcat('binSpkCountSTRCTX',p.Results.fileInfo);
save(fullfile(p.Results.filePath,saveName),'-struct','binSpkCountSTRCTX') % save the fields of the structure separately 
save(fullfile(p.Results.filePath,saveName), 'evt', 'jkvt', '-append') % append the behavioral timestamps
clearvars saveName binSpkCountSTRCTX pullStarts pullStops trStart trEnd reward stmLaser tagLaser rStartToPull rStopToPull 

%% Individual unit raster plot
% cUnits = find(cell2mat(pullStarts.isStr)==0); 
% sUnits = find(cell2mat(pullStarts.isStr)==1); 
% %sUnitN = 1; 
% spikeRasterGramm( [3e3 2e3], {'pullStarts'}, [2e3 2e3], binSpkCount.pullStarts.SpkTimes{sUnits(sUnitN)});
% spikeRasterGramm( [3e3 2e3], {'rStartToPull'}, [2e3 2e3], binSpkCount.rStartToPull.SpkTimes{sUnits(sUnitN)});
% sUnitN = sUnitN+1; 
%print( fullfile(filePath,'Figure',strcat(fileInfo,'_',sprintf('unit#%d',unitNumb),'pullStart')), '-dpdf','-painters', '-bestfit')

%unit = 139; % str unit 128 (laser activated)
%spikeRasterGramm( [1e3 3e3], {'reach','stim','stimReach'}, binSpkCountSTRReach(unit).SpkTimes, binSpkCountSTRstmLaser(unit).SpkTimes,binSpkCountSTRstmReach(unit).SpkTimes );


%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
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







