function PSTH_rasters_reachDecodeSTR( filePath, fileName, fileInfo )
%PSTH_rasters takes behavioral timestamps stored in BehVariables.mat and
% generates psths aligned to each event and saves the outcome psths in the
% filePath. To run PSTH_rasters 'behaviorTimestamps.m' must be run first
% and its outcome 'BehVariables.mat' must exist in the filePath. 
% Modified on 6/18/18 to save the combined (all imec channels)
% 'binSpkCountStrCtx*.mat'
% Modified in Jan/19 to change the baseline period for tagLaser PSTHs from
% reachStart to tagLaser. 

%% Load files 
load(fullfile(filePath,fileName),'spkTimesCellSTR') % load spkTimesCellCtx
load(fullfile(filePath,'BehVariablesReachDecode.mat'),'ts') % load 

% binned spike count STR 
reachMore = psthBINcell( fileInfo, 'STR', spkTimesCellSTR, ts.reachStartMore', ts.reachStartMore', 1, [2e3 3e3], -1, false );    % entire reachStart 
noReach = psthBINcell( fileInfo, 'STR', spkTimesCellSTR, ts.noReachTime', ts.noReachTime', 1, [2e3 3e3], -1, false );    % entire reachStart 

binSpkCountSTRreachDecode.reachMore = reachMore; 
binSpkCountSTRreachDecode.noReach = noReach; 

saveNameSTR = strcat('binSpkCountSTRreachDecode',fileInfo);
save(saveNameSTR, '-struct', 'binSpkCountSTRreachDecode') % save the fields of the structure separately 

end


