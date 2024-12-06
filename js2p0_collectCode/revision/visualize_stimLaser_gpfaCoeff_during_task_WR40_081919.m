% For visualization, this needs to be run first: js2p0_tbytSpkHandJsPreprocess_50ms_stim_parse(filePath)
filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles'; 
load(fullfile(filePath, strcat('js2p0_tbytSpkHandJsTrjBin_50ms_stimParse_', 'WR40_081919')),'ss', 'jkvt'); 
%load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI') 
filePath_gpfa = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/gpfa_ctx'; 
load(fullfile(filePath_gpfa, 'gpfaRezCtx.mat'), 'gpfaRezCtx'); 

figSaveDir = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/Figure'; 

%%
ctx_base_mean = nanmean(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBCtx}, 'UniformOutput', false))), 2); 
ctx_base_std = nanstd(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBCtx}, 'UniformOutput', false))), 0, 2); 

for t = 1:length(ss)
    ss(t).utbCtxZ = (full(ss(t).unitTimeBCtx)-repmat(ctx_base_mean, 1, size(ss(t).unitTimeBCtx, 2)))...
        ./repmat(ctx_base_std, 1, size(ss(t).unitTimeBCtx, 2));
end

%% identify cells with significant silencing 
x = -975:50:3975; 
pre_stim = [];
post_stim = []; 
post_stim_min = []; 
breakthrough = []; 
reach = []; 
for t = 1:length(ss)
    if ~isempty(ss(t).utbCtxStimAlignZ) && (isempty(ss(t).rStartRelToStim) || ss(t).rStartRelToStim>2000)
       pre_stim = [pre_stim, nanmean(ss(t).utbCtxStimAlignZ(:, 1:20), 2)]; 
       post_stim_min = [post_stim, nanmin(ss(t).utbCtxStimAlignZ(:, 21:40), [], 2)]; 
       post_stim = [post_stim, nanmean(ss(t).utbCtxStimAlignZ(:, 21:40), 2)]; 
    end

    if ~isempty(ss(t).utbCtxStimAlignZ) && ~isempty(ss(t).rStartRelToStim) && ss(t).rStartRelToStim<3500
        breakthrough = [breakthrough, nanmean(ss(t).utbCtxStimAlignZ(:, x>ss(t).rStartRelToStim), 2)]; 
    end

    if strcmpi(ss(t).trialType, 'sp') && isempty(ss(t).tLaserStart)
        reach = [reach, nanmean(full(ss(t).utbCtxZ), 21:40)]; 
    end
end

for i = 1:size(pre_stim, 1)
    ttestRez(i, 1) = ttest(pre_stim(i,:), post_stim(i,:), 'Alpha', 0.01); 
    if nanmean(pre_stim(i,:))<nanmean(post_stim(i,:)) || isnan(ttestRez(i, 1))
        ttestRez(i, 1) = 0;
    end
end

rez.magInh = nanmean(post_stim_min, 2)-nanmean(pre_stim, 2); 
rez.sigInhI = logical(ttestRez); 
rez.magBT = nanmean(breakthrough, 2)-nanmean(pre_stim, 2); 
rez.magReach = nanmean(reach, 2)-nanmean(pre_stim, 2); 
rez.gpfaCoeff = gpfaRezCtx.faParams.L; 

%scatter_row_by_row_fillLogic([rez.magInh, rez.magBT, rez.magReach], [-1 3],  logical(ttestRez))
save(fullfile(filePath, 'magnitude_inhibition_breakThrough_gpfa'), 'rez')






