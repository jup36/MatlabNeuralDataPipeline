%This script loads 1) behavioral timestamps (generated by
% BehaviorTimestamps*.m), 2) jrc.mat (jrclust), 
% 3) meta file containing the recording info (spikeGLX) to generate PSTHs.
% Make sure to have 1)*ap.bin, 2)*ap.meta, 3)BehVariables.mat, 4)*_jrc.mat in the working path. 

clear all; clear functions; clc

%% Raw voltage traces from bin and meta files
addpath(genpath('/Volumes/RAID2/parkj/MATLAB')) % folder with relevant matlab codes
addpath(genpath('/Users/parkj/Desktop/SpikeGLX-master 2/MATLAB-SDK')) % the folder with spikeGLX matlab scripts by Bill Karsh 
addpath(genpath('/Volumes/RAID2/parkj/jrclust_repo/JRCLUST')) % add jrclust-repo

%% Load files 
path    = '/Volumes/RAID2/parkj/NeuralData/IT03_Ldms_M1_013118/Matfiles';    % file directory for matfiles 
cd(path)   % change directory to the current data folder

geometry = getimec3opt3geom; % get the geometry of the imec3opt3 probe
meta  = getmeta;             % get meta file using the helper function getmeta
[S_clu,viTime_spk,viSite_spk] = getjrcmatVar; % get the structure variable containing cluster info
load('BehVariables', 'ts')   % load behavioral timestamps (reach start, stop, and position/velocity data for each reach)

%% File information
fileInfo       = 'IT03_013118';   % to put the essential file information (animal ID and date of data collection) into the structure for each single units. 
dvCosConvert   = cos(10/180*pi);  % if probe was angled, probe coordinates need to be corrected 
probeDepth     = [3985];          % ventral coordinate (in um) of the dms probe relative to the brain surface 
strCtx         = true;            % boolean to indicate recording from both striatum and cortex on the same axis (true, e.g. imec probe) 
strCtxBorder   = 2000;            % dv relative to the brain surface to be counted as a striatal unit (striatal units should have a dv greater than this)
numbBehSite    = 64;              % the number of NI board sites used for behavioral data collection (default=64)
numbSiteProbe  = [384];           % the number of sites per probe
numbProbe      = unique([length(probeDepth) length(numbSiteProbe)]); % the number of probes used

if length(probeDepth)~=length(numbSiteProbe)
    error('Check the variables; probeDepth & numbChProbe!')
end

% Assign probe id to channels (as multiple probes might be used)
whichProbe = zeros(sum(numbSiteProbe),1); % probe identity for each channel  

probeId = 1; % one-based (first probe gets 1)
for s = 1:sum(numbSiteProbe) % increment all probe channels
    if s <= sum(numbSiteProbe(1:probeId))
        whichProbe(s,1) = probeId;
    elseif s > sum(numbSiteProbe(1:probeId))
        whichProbe(s,1) = probeId+1;
        probeId = probeId+1;
    end
end
clearvars probeId s 

%% Process spike times data and generate PSTH and Rasters
% get spike times from the csv file
spkTimes = struct; % the structure to contain spike times
for u = 1:S_clu.nClu % increment valid clusters (units)
    clusId = ['c' num2str(u)];    % cluster name string e.g. c1, c2, c3 ...
    
    spkIdx = S_clu.cviSpk_clu{u}; % spikeIDs of the current cluster 
    spkTimes(u).spkTimes = double(viTime_spk(spkIdx))/30; % divided by the sampling rate: 30kHz
    spkTimes(u).clusId   = u; % cluster ID
    spkTimes(u).maxSite  = mode(double(viSite_spk(spkIdx))); % assign the current cluster to a site of the probe
    spkTimes(u).geometry(1) = geometry(spkTimes(u).maxSite,1);        % horizontal geometry
    spkTimes(u).geometry(2) = probeDepth(whichProbe(spkTimes(u).maxSite))-geometry(spkTimes(u).maxSite,2); % vertical geometry
    spkTimes(u).geometry(2) = spkTimes(u).geometry(2).*dvCosConvert;  % corrected dv (original dv multiplied by cosine(probe angle))
    
    if strCtx && strcmp(meta.typeThis, 'imec') % in case recording from both str and ctx using imec probe
        spkTimes(u).isStr = spkTimes(u).geometry(2) >= strCtxBorder;  % logical for striatum 
    end
    
end
clearvars u 

spkTimesCell    = struct2cell(spkTimes'); % the entire spike times converted into a cell 
if strcmp(meta.typeThis, 'imec')          % in case recording via imec 
    spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(5,:))==0); % the CTX spike times cell (1st probe)
    spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(5,:))==1); % the STR spike times cell (2nd probe)
elseif strcmp(meta.typeThis, 'nidq')      % in case recording from one or two apig probes
    spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(3,:))<=64); % the CTX spike times cell (1st probe)
    spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(3,:))>64);  % the STR spike times cell (2nd probe)
end

psthPlotFlag   = false; % boolean - to plot psth or not

% binned spike count CTX reachStart 
binSpkCountCTXReach    = psthBIN( fileInfo, 'M1', spkTimesCellCTX, ts.reachStart', ts.reachStart', 1, [1e3 3e3], -1, psthPlotFlag );    % entire reachStart 
binSpkCountCTXrwd      = psthBIN( fileInfo, 'M1', spkTimesCellCTX, ts.reward', ts.reachStart', 1, [3e3 1e3], -1, psthPlotFlag );        % entire rewardDelivery
binSpkCountCTXstmLaser = psthBIN( fileInfo, 'M1', spkTimesCellCTX, ts.stmLaser', ts.reachStart', 1, [1e3 3e3], -1, psthPlotFlag );      % laser stim trials
binSpkCountCTXtagLaser = psthBIN( fileInfo, 'M1', spkTimesCellCTX, ts.tagLaser', ts.reachStart', 1, [5e3 5e3], -1, psthPlotFlag );      % laser tag trials
binSpkCountCTXstmReach = psthBIN( fileInfo, 'M1', spkTimesCellCTX, ts.stmReachStart', ts.reachStart', 1, [1e3 3e3], -1, psthPlotFlag ); % reachStart with stimulation (completed reaches even during stimulation)

% binned spike count STR reachStart
binSpkCountSTRReach    = psthBIN( fileInfo, 'DMS', spkTimesCellSTR, ts.reachStart', ts.reachStart', 1, [1e3 3e3], -1, psthPlotFlag );    % entire reachStart 
binSpkCountSTRrwd      = psthBIN( fileInfo, 'DMS', spkTimesCellSTR, ts.reward', ts.reachStart', 1, [3e3 1e3], -1, psthPlotFlag );        % entire rewardDelivery
binSpkCountSTRstmLaser = psthBIN( fileInfo, 'DMS', spkTimesCellSTR, ts.stmLaser', ts.reachStart', 1, [1e3 3e3], -1, psthPlotFlag );      % laser stim trials
binSpkCountSTRtagLaser = psthBIN( fileInfo, 'DMS', spkTimesCellSTR, ts.tagLaser', ts.reachStart', 1, [5e3 5e3], -1, psthPlotFlag );      % laser tag trials
binSpkCountSTRstmReach = psthBIN( fileInfo, 'DMS', spkTimesCellSTR, ts.stmReachStart', ts.reachStart', 1, [1e3 3e3], -1, psthPlotFlag ); % reachStart with stimulation (completed reaches even during stimulation)

% save files
saveNameCTX = strcat('binSpkCountCTX',fileInfo);
saveNameSTR = strcat('binSpkCountSTR',fileInfo);
save(saveNameCTX, 'binSpkCountCTXReach', 'binSpkCountCTXrwd', 'binSpkCountCTXstmLaser', 'binSpkCountCTXtagLaser', 'binSpkCountCTXstmReach', 'spkTimesCellCTX', 'meta') % save the data
save(saveNameSTR, 'binSpkCountSTRReach', 'binSpkCountSTRrwd', 'binSpkCountSTRstmLaser', 'binSpkCountSTRtagLaser', 'binSpkCountSTRstmReach', 'spkTimesCellSTR', 'meta') % save the data

%% Individual unit raster plot
spikeRasterGramm( [1e3 3e3], {'reach','stimReach'}, binSpkCountCTXReach(11).SpkTimes, binSpkCountCTXstmReach(11).SpkTimes );

unit = 139; % str unit 128 (laser activated)
spikeRasterGramm( [1e3 3e3], {'reach','stim','stimReach'}, binSpkCountSTRReach(unit).SpkTimes, binSpkCountSTRstmLaser(unit).SpkTimes,binSpkCountSTRstmReach(unit).SpkTimes );




