function PSTH_rasters_js2p0_trJsReady( filePath )
%'PSTH_rasters_js2p0_trJsReady.m' gets PSTHs around the trial joystick ready
% events and append the data to the existing 'binSpkCount*' files.

%% Ctx
bsc_ctx = dir(fullfile(filePath, '**/*binSpkCountCTXWR*'));
if ~isempty(bsc_ctx)
    load(fullfile(bsc_ctx(1).folder, bsc_ctx(1).name), 'spkTimesCell', 'jkvt', 'p')
    [trJsReady_indices] = trJsReady_typify(jkvt);
    trJsReady = psthBINcell( p.Results.fileInfo, 'M1', spkTimesCell, [jkvt.trJsReady]', [jkvt.trJsReady]'-1000, 1, [2e3 3e3], -1, p.Results.psthPlotFlag );
    trJsReady.rStartToPullIdx = trJsReady_indices.rStartToPullIdx; 
    trJsReady.pullStartsIdx = trJsReady_indices.pullStartsIdx; 
    trJsReady.rewardIdx = trJsReady_indices.rewardIdx; 
    save(fullfile(bsc_ctx(1).folder, bsc_ctx(1).name), 'trJsReady', '-append')
    clearvars -except filePath jkvt
end

%% Str
bsc_str = dir(fullfile(filePath, '**/*binSpkCountSTRWR*'));
if ~isempty(bsc_str)
    load(fullfile(bsc_str(1).folder, bsc_str(1).name), 'spkTimesCell', 'p')
    
    if ~exist('jkvt', 'var')
        load(fullfile(bsc_str(1).folder, bsc_str(1).name), 'spkTimesCell', 'jkvt', 'p')
    end
    
    [trJsReady_indices] = trJsReady_typify(jkvt);
    trJsReady = psthBINcell( p.Results.fileInfo, 'DMS', spkTimesCell, [jkvt.trJsReady]', [jkvt.trJsReady]'-1000, 1, [2e3 3e3], -1, p.Results.psthPlotFlag );
    trJsReady.rStartToPullIdx = trJsReady_indices.rStartToPullIdx; 
    trJsReady.pullStartsIdx = trJsReady_indices.pullStartsIdx; 
    trJsReady.rewardIdx = trJsReady_indices.rewardIdx; 
    save(fullfile(bsc_str(1).folder, bsc_str(1).name), 'trJsReady', '-append')
    clearvars -except filePath jkvt
end

%% CtxStr
bsc_ctx_str = dir(fullfile(filePath, '**/*binSpkCountSTRCTXWR*'));
if ~isempty(bsc_ctx_str)
    load(fullfile(bsc_ctx_str(1).folder, bsc_ctx_str(1).name), 'spkTimesCell', 'p')
    
    if ~exist('jkvt', 'var')
        load(fullfile(bsc_ctx_str(1).folder, bsc_ctx_str(1).name), 'spkTimesCell', 'jkvt', 'p')
    end
    
    [trJsReady_indices] = trJsReady_typify(jkvt);
    trJsReady = psthBINcell( p.Results.fileInfo, 'DMSM1', spkTimesCell, [jkvt.trJsReady]', [jkvt.trJsReady]'-1000, 1, [2e3 3e3], -1, p.Results.psthPlotFlag );
    trJsReady.rStartToPullIdx = trJsReady_indices.rStartToPullIdx; 
    trJsReady.pullStartsIdx = trJsReady_indices.pullStartsIdx; 
    trJsReady.rewardIdx = trJsReady_indices.rewardIdx; 
    save(fullfile(bsc_ctx_str(1).folder, bsc_ctx_str(1).name), 'trJsReady', '-append')
    clearvars -except filePath jkvt
end

%% Cg
bsc_cg = dir(fullfile(filePath, '**/*binSpkCountCgWR*'));
if ~isempty(bsc_cg)
    load(fullfile(bsc_cg(1).folder, bsc_cg(1).name), 'spkTimesCell', 'p')
    
    if ~exist('jkvt', 'var')
        load(fullfile(bsc_cg(1).folder, bsc_cg(1).name), 'spkTimesCell', 'jkvt', 'p')
    end
    
    [trJsReady_indices] = trJsReady_typify(jkvt);
    trJsReady = psthBINcell( p.Results.fileInfo, 'Cg', spkTimesCell, [jkvt.trJsReady]', [jkvt.trJsReady]'-1000, 1, [2e3 3e3], -1, p.Results.psthPlotFlag );
    trJsReady.rStartToPullIdx = trJsReady_indices.rStartToPullIdx; 
    trJsReady.pullStartsIdx = trJsReady_indices.pullStartsIdx; 
    trJsReady.rewardIdx = trJsReady_indices.rewardIdx; 
    save(fullfile(bsc_cg(1).folder, bsc_cg(1).name), 'trJsReady', '-append')
    clearvars -except filePath jkvt
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
function [trJsReady] = trJsReady_typify(jkvt)
trJsReady_trI = find(~cellfun(@isempty, {jkvt.trJsReady}));
rStartToPull_trI = find(~cellfun(@isempty, {jkvt.rStartToPull}));
pullStarts_trI = find(~cellfun(@isempty, {jkvt.pullStarts}));
reward_trI = find(cellfun(@(a) a == 1, {jkvt.rewarded}));

trJsReady.rStartToPullIdx = ismember(trJsReady_trI, rStartToPull_trI);
trJsReady.pullStartsIdx = ismember(trJsReady_trI, pullStarts_trI);
trJsReady.rewardIdx = ismember(trJsReady_trI, reward_trI);
end



end









