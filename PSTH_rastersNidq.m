function PSTH_rastersNidq( filePath, fileInfo, probeDepth, varargin )
%PSTH_rasters takes behavioral timestamps stored in BehVariables.mat and
% generates psths aligned to each event and saves the outcome psths in the
% filePath. To run PSTH_rasters 'behaviorTimestamps.m' must be run first
% and its outcome 'BehVariables.mat' must exist in the filePath. 

p = parse_input_psthNidq(filePath, fileInfo, probeDepth, varargin ); % parse input
%p = parse_input_psthNidq(filePath,'WR23_091517',[1724 2930],{}) % when running line-by-line

%% Load files 
cd(p.Results.filePath)   % change directory to the data folder

behFile = dir(fullfile(p.Results.filePath,'BehVariables.mat')); % look for 'BehVariables.mat' file    

if length(behFile)>1 || isempty(behFile)    
    error('File could not be found or multiple BehVariables.mat files exist!');
else
    load(behFile.name,'ts') % load behavioral timestamps 
end

geometry = [];
for i = 1:length(p.Results.numbSiteProbe) % increment the number of probes
    tempgeometry = getHH2AcuteGeom; % get the geometry of the imec3opt3 probe
    geometry = [geometry; tempgeometry ];
end

meta  = getmeta;            % get meta file using the helper function getmeta
[S_clu,viTime_spk,viSite_spk] = getjrcmatVar; % get the structure variable containing cluster info

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
spkTimes = struct; % the structure to contain spike times
for u = 1:S_clu.nClu % increment valid clusters (units)
    
    spkIdx = S_clu.cviSpk_clu{u}; % spikeIDs of the current cluster 
    
    if strcmp(meta.typeThis, 'imec')
        spkTimes(u).spkTimes = double(viTime_spk(spkIdx))/(str2double(meta.imSampRate)/1000); % divided by the sampling rate: 30kHz (str2num(meta.imSampRate)/1000)
    elseif strcmp(meta.typeThis, 'nidq')
        spkTimes(u).spkTimes = double(viTime_spk(spkIdx))/(str2double(meta.niSampRate)/1000); % divided by the sampling rate: 25kHz (str2num(meta.niSampRate)/1000)
    end
    
    spkTimes(u).clusId   = u; % cluster ID
    spkTimes(u).maxSite  = mode(double(viSite_spk(spkIdx))); % assign the current cluster to a site of the probe
    spkTimes(u).geometry(1) = geometry(spkTimes(u).maxSite,1);        % horizontal geometry
    spkTimes(u).geometry(2) = p.Results.probeDepth(whichProbe(spkTimes(u).maxSite))-geometry(spkTimes(u).maxSite,2); % vertical geometry
    spkTimes(u).geometry(2) = spkTimes(u).geometry(2).*dvCosConvert;  % corrected dv (original dv multiplied by cosine(probe angle))
    
    if strcmp(meta.typeThis, 'imec') % in case recording from both str and ctx using imec probe
       spkTimes(u).isStr = spkTimes(u).geometry(2) >= p.Results.strCtxBorder;  % logical for striatum 
    elseif strcmp(meta.typeThis, 'nidq')
       spkTimes(u).isStr = [];  %
    end
    
    spkTimes(u).meanWF = S_clu.tmrWav_raw_clu(:,S_clu.viSite_clu(u),u); % mean raw waveforms for each cluster
end
clearvars u 

spkTimesCell    = struct2cell(spkTimes'); % the entire spike times converted into a cell 
if strcmp(meta.typeThis, 'imec')          % in case recording via imec  
    spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(5,:))==0); % the CTX spike times cell (cells from sites above the strCtxBorder)
    spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(5,:))==1); % the STR spike times cell (cells from sites below the strCtxBorder )
elseif strcmp(meta.typeThis, 'nidq')      % in case recording from one or two apig probes
    spkTimesCellCTX = spkTimesCell(:,cell2mat(spkTimesCell(3,:))<=64); % the CTX spike times cell (1st probe)
    spkTimesCellSTR = spkTimesCell(:,cell2mat(spkTimesCell(3,:))>64);  % the STR spike times cell (2nd probe)
end

% binned spike count CTX reachStart 
binSpkCountCTX.reach    = psthBINcell( p.Results.fileInfo, 'Ctx', spkTimesCellCTX, ts.reachStart', ts.reachStart', 1, p.Results.reachWin, -1, p.Results.psthPlotFlag );    % entire reachStart 
binSpkCountCTX.reward   = psthBINcell( p.Results.fileInfo, 'Ctx', spkTimesCellCTX, ts.reward', ts.reachStart', 1, p.Results.rewardWin, -1, p.Results.psthPlotFlag );       % entire rewardDelivery
if isfield(ts,'stmLaser')
    binSpkCountCTX.stmLaser = psthBINcell( p.Results.fileInfo, 'Ctx', spkTimesCellCTX, ts.stmLaser', ts.reachStart', 1, p.Results.reachWin, -1, p.Results.psthPlotFlag );      % laser stim trials
end

if isfield(ts,'tagLaser')
    binSpkCountCTX.tagLaser = psthBINcell( p.Results.fileInfo, 'Ctx', spkTimesCellCTX, ts.tagLaser', ts.reachStart', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag );   % laser tag trials
end

if isfield(ts,'stmReachStart')
    binSpkCountCTX.stmReach = psthBINcell( p.Results.fileInfo, 'Ctx', spkTimesCellCTX, ts.stmReachStart', ts.reachStart', 1, p.Results.reachWin, -1, p.Results.psthPlotFlag ); % reachStart with stimulation (completed reaches even during stimulation)
end
    
% binned spike count STR reachStart
binSpkCountSTR.reach    = psthBINcell( p.Results.fileInfo, 'Str', spkTimesCellSTR, ts.reachStart', ts.reachStart', 1, p.Results.reachWin, -1, p.Results.psthPlotFlag );    % entire reachStart 
binSpkCountSTR.reward   = psthBINcell( p.Results.fileInfo, 'Str', spkTimesCellSTR, ts.reward', ts.reachStart', 1, p.Results.rewardWin, -1, p.Results.psthPlotFlag );       % entire rewardDelivery

if isfield(ts,'stmLaser')
    binSpkCountSTR.stmLaser = psthBINcell( p.Results.fileInfo, 'Str', spkTimesCellSTR, ts.stmLaser', ts.reachStart', 1, p.Results.reachWin, -1, p.Results.psthPlotFlag );      % laser stim trials
end

if isfield(ts,'stmReachStart')
    binSpkCountSTR.tagLaser = psthBINcell( p.Results.fileInfo, 'Str', spkTimesCellSTR, ts.tagLaser', ts.reachStart', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag );   % laser tag trials
end

if isfield(ts,'stmReachStart')
    binSpkCountSTR.stmReach = psthBINcell( p.Results.fileInfo, 'Str', spkTimesCellSTR, ts.stmReachStart', ts.reachStart', 1, p.Results.reachWin, -1, p.Results.psthPlotFlag ); % reachStart with stimulation (completed reaches even during stimulation)
end

% to save meta & p info
binSpkCountCTX.meta = meta; 
binSpkCountCTX.p = p; 
binSpkCountSTR.meta = meta; 
binSpkCountSTR.p = p; 
binSpkCountCTX.spkTimesCellCTX = spkTimesCellCTX; % just to save this
binSpkCountSTR.spkTimesCellSTR = spkTimesCellSTR; % just to save this

% save files
saveNameCTX = strcat('binSpkCountCTX',p.Results.fileInfo);
saveNameSTR = strcat('binSpkCountSTR',p.Results.fileInfo);

save(saveNameCTX, '-struct', 'binSpkCountCTX') % save the fields of the structure separately 
save(saveNameCTX, 'ts', '-append') % append the behavioral timestamps
save(saveNameSTR, '-struct', 'binSpkCountSTR') % save the fields of the structure separately
save(saveNameSTR, 'ts', '-append') % append the behavioral timestamps
 
%% Individual unit raster plot
%spikeRasterGramm( [1e3 3e3], {'reach','stimReach'}, binSpkCountCTXReach(11).SpkTimes, binSpkCountCTXstmReach(11).SpkTimes );

%unit = 139; % str unit 128 (laser activated)
%spikeRasterGramm( [1e3 3e3], {'reach','stim','stimReach'}, binSpkCountSTRReach(unit).SpkTimes, binSpkCountSTRstmLaser(unit).SpkTimes,binSpkCountSTRstmReach(unit).SpkTimes );


%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_psthNidq( filePath, fileInfo, probeDepth, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function 
        %parse input, and extract name-value pairs for the main function
        % 'PSTH_rasters'
                
        default_probeAngle = 10;       % default probe angle in degree
        default_numbSiteProbe = [64];  % default number of sites per probe(s) (can be an arrary, if there are multiple probes used
        default_psthPlotFlag = false;  % default logic indicating psth draw or not 
        default_reachWin = [2e3 2e3];  % default time window for reach psth
        default_rewardWin = [3e3 1e3]; % default time window for reward psth
        default_tagLaserWin = [5e3 5e3]; % default time window for tagLaser psth
        default_area1 = 'Ctx';
        default_area2 = 'Str';
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addRequired(p,'fileInfo');
        addRequired(p,'probeDepth'); % multiple inputs possible when multiple probes were used

        addParameter(p,'probeAngle', default_probeAngle)        
        addParameter(p,'numbSiteProbe', default_numbSiteProbe)
        addParameter(p,'psthPlotFlag', default_psthPlotFlag)
        addParameter(p,'reachWin', default_reachWin)
        addParameter(p,'rewardWin', default_rewardWin)
        addParameter(p,'tagLaserWin', default_tagLaserWin)
        addParameter(p,'area1', default_area1)
        addParameter(p,'area2', default_area2)
        
        parse(p,filePath, fileInfo, probeDepth, vargs{:})
        
    end


end


%%%%%%%%%%%%%%%%%%%%%%%
%  NON-NESTED FUCNTIONS
%%%%%%%%%%%%%%%%%%%%%%%
function [geometry] = getHH2AcuteGeom
    vi32 = 1:32;
    geometry = zeros(64, 2);
    geometry(vi32,2) = 25*(vi32-1);
    geometry(vi32+32,1) = 250;
    geometry(vi32+32,2) = geometry(vi32,2);
end


function [ geometry ] = getimec3opt3geom
%This function returns the geometry of the imec3 option3 probe. 

% Order of the probe sites in the recording file
channels = 1:384;

% Site location in micrometers (x and y)
geometry = zeros(384, 2);
viHalf = 0:(384/2-1);
geometry(1:2:end,2) = viHalf * 20;
geometry(2:2:end,2) = geometry(1:2:end,2);
geometry(1:4:end,1) = 16; %0 before
geometry(2:4:end,1) = 48; %32 before
geometry(3:4:end,1) = 0;  %16 before
geometry(4:4:end,1) = 32; %48 before

% Reference sites are being excluded
ref_sites = [37 76 113 152 189 228 265 304 341 380];
channels(ref_sites) = []; 
geometry(ref_sites,:) = []; 

% Recording contact pad size in micrometers. Height x width
pad = [12 12];

% Default prm
um_per_pix = 20;
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







