function PSTH_rasters_js2p0_rec( filePath, fileInfo, probeDepth, varargin )
%PSTH_rasters_js2p0_rec replaces 'PSTH_rasters_js2p0.m' to generate
% binSpkCount*.mat files that are lost due to the hard disk (beefcake)
% failure. Main inputs are 'evtIndices.mat',
% 'js2p0_tbytSpkHandJsTrjBin_*.mat', which has 'jkvt' and 'spikeTimeCells'
% that are intact. 
% Modified on Apr/21

p = parse_input_psth_js2p0(filePath, fileInfo, probeDepth, varargin); % parse input
%p = parse_input_psth_js2p0('/Volumes/8TB/Junchol_Data/JS2p0/WR40_081919', 'WR40_081919_rec', 4172, {'probeAngle', 10, 'numbSiteProbe',384, 'laserUsed', true}) % when running line-by-line

% get meta 
disp('Select the meta file!!')
[metaFileSel,metaPathSel] = uigetfile('*.meta',p.Results.filePath); 
meta = ReadMeta(metaFileSel, metaPathSel); % read out the meta file

%% Load files 
cd(p.Results.filePath)   % change directory to the data folder

% get behavioral data
evtFile = dirsub(p.Results.filePath, 'evtIndices.mat'); 
if length(evtFile)>1 || isempty(evtFile)
    disp('Select evtIndices.mat!')
    [evtFileSelect,evtPathSelect] = uigetfile('*.mat',p.Results.filePath);
    evt = load(fullfile(evtPathSelect,evtFileSelect),'evtIdx1k');
    evt = evt.('evtIdx1k'); 
    clearvars evtIdx1k
else
    evt = load(fullfile(evtFile{1}.folder,evtFile{1}.name),'evtIdx1k');
    evt = evt.('evtIdx1k'); 
    clearvars evtIdx1k
end

jkvtStcFile = dirsub(p.Results.filePath, 'js2p0_tbytSpkHandJsTrjBin_*.mat'); 
if length(jkvtStcFile)>1 || isempty(jkvtStcFile)
    disp('js2p0_tbytSpkHandJsTrjBin_*.mat!')
    [jkvtStcFileSelect,jkvtStcPathSelect] = uigetfile('*.mat',p.Results.filePath);
    load(fullfile(jkvtStcPathSelect,jkvtStcFileSelect),'jkvt');
    load(fullfile(jkvtStcPathSelect,jkvtStcFileSelect),'spkTimesCell'); 
else
    load(fullfile(jkvtStcFile{1}.folder,jkvtStcFile{1}.name),'jkvt');
    load(fullfile(jkvtStcFile{1}.folder,jkvtStcFile{1}.name),'spkTimesCell'); 
end

isStr = [spkTimesCell{5,:}]; 
spkTimesCellStrCtx = spkTimesCell; 
spkTimesCellSTR = spkTimesCellStrCtx(:,isStr); 
spkTimesCellCTX = spkTimesCellStrCtx(:,~isStr); 

%% Specify all the time points to align neural data onto 
% define pull Starts/Stops
evt.trJsReady = [jkvt.trJsReady]'; % joystick ready timePoint per trial
pullTrsIdx = cellfun(@(c)strcmpi(c,'sp'), {jkvt.trialType}); % successfull pull trials

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
if p.Results.laserUsed
    stmLaser       = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, evt.stimLaserRiseIdx', evt.stimLaserRiseIdx', 1, [1e3 5e3], -1, p.Results.psthPlotFlag ); % laser stim trials
    tagLaser       = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCellCTX, evt.tagLaserRiseIdx', evt.tagLaserRiseIdx', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag ); % laser tag trials
    
    binSpkCountCTX.stmLaser = stmLaser;
    binSpkCountCTX.tagLaser = tagLaser;
end
binSpkCountCTX.pullStarts = pullStarts; 
binSpkCountCTX.pullStops  = pullStops; 
binSpkCountCTX.rStartToPull = rStartToPull; 
binSpkCountCTX.rStopToPull = rStopToPull; 
binSpkCountCTX.reward = reward; 
binSpkCountCTX.meta = meta; 
binSpkCountCTX.p = p; 
binSpkCountCTX.spkTimesCell = spkTimesCellCTX; % just to save the cell

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

if p.Results.laserUsed
    stmLaser       = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, evt.stimLaserRiseIdx', evt.stimLaserRiseIdx', 1, [1e3 5e3], -1, p.Results.psthPlotFlag ); % laser stim trials
    tagLaser       = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCellSTR, evt.tagLaserRiseIdx', evt.tagLaserRiseIdx', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag ); % laser tag trials

    binSpkCountSTR.stmLaser = stmLaser; 
    binSpkCountSTR.tagLaser = tagLaser; 
end

binSpkCountSTR.pullStarts = pullStarts; 
binSpkCountSTR.pullStops  = pullStops; 
binSpkCountSTR.rStartToPull = rStartToPull; 
binSpkCountSTR.rStopToPull = rStopToPull; 
binSpkCountSTR.reward = reward; 
binSpkCountSTR.meta = meta; 
binSpkCountSTR.p = p; 
binSpkCountSTR.spkTimesCell = spkTimesCellSTR; % just to save the cell

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
if p.Results.laserUsed
    stmLaser       = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, evt.stimLaserRiseIdx', evt.stimLaserRiseIdx', 1, [1e3 5e3], -1, p.Results.psthPlotFlag ); % laser stim trials
    tagLaser       = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCellStrCtx, evt.tagLaserRiseIdx', evt.tagLaserRiseIdx', 1, p.Results.tagLaserWin, -1, p.Results.psthPlotFlag ); % laser tag trials
    
    binSpkCountSTRCTX.stmLaser = stmLaser;
    binSpkCountSTRCTX.tagLaser = tagLaser;
end
binSpkCountSTRCTX.pullStarts = pullStarts; 
binSpkCountSTRCTX.pullStops  = pullStops; 
binSpkCountSTRCTX.rStartToPull = rStartToPull; 
binSpkCountSTRCTX.rStopToPull = rStopToPull; 
binSpkCountSTRCTX.reward = reward; 
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

function subdirinfo = dirsub(filePath, fileName)

dirinfo = dir(filePath);
dirinfo(~[dirinfo.isdir]) = [];

subdirinfo = {};
for K = 1 : length(dirinfo)
    thisdir = fullfile(dirinfo(K).folder,dirinfo(K).name);
    subdirinfo = [subdirinfo; dir(fullfile(thisdir, fileName))]; 
end

end







