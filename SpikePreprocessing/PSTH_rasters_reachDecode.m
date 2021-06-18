function PSTH_rasters_reachDecode( filePath, fileName, fileInfo )
%PSTH_rasters takes behavioral timestamps stored in BehVariables.mat and
% generates psths aligned to each event and saves the outcome psths in the
% filePath. To run PSTH_rasters 'behaviorTimestamps.m' must be run first
% and its outcome 'BehVariables.mat' must exist in the filePath. 
% Modified on 6/18/18 to save the combined (all imec channels)
% 'binSpkCountStrCtx*.mat'
% Modified in Jan/19 to change the baseline period for tagLaser PSTHs from
% reachStart to tagLaser. 

%% Load files 
load(fullfile(filePath,fileName),'spkTimesCellCTX') % load spkTimesCellCtx
load(fullfile(filePath,'BehVariablesReachDecode.mat'),'ts') % load 

% binned spike count CTX 
reachMore = psthBINcell( fileInfo, 'M1', spkTimesCellCTX, ts.reachStartMore', ts.reachStartMore', 1, [2e3 3e3], -1, false );    % entire reachStart 
noReach = psthBINcell( fileInfo, 'M1', spkTimesCellCTX, ts.noReachTime', ts.noReachTime', 1, [2e3 3e3], -1, false );    % entire reachStart 

binSpkCountCTXreachDecode.reachMore = reachMore; 
binSpkCountCTXreachDecode.noReach = noReach; 

saveNameCTX = strcat('binSpkCountCTXreachDecode',fileInfo);
save(saveNameCTX, '-struct', 'binSpkCountCTXreachDecode') % save the fields of the structure separately 

end


