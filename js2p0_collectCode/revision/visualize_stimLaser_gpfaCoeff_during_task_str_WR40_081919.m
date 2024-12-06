% For visualization, this needs to be run first: js2p0_tbytSpkHandJsPreprocess_50ms_stim_parse(filePath)
filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles'; 
load(fullfile(filePath, strcat('js2p0_tbytSpkHandJsTrjBin_50ms_stimParse_', 'WR40_081919')),'ss', 'jkvt'); 
%load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI') 
%filePath_gpfa = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/gpfa_str'; 
%load(fullfile(filePath_gpfa, 'gpfaRezStr.mat'), 'gpfaRezStr'); 

figSaveDir = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/Figure'; 

%%
str_base_mean = nanmean(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBStr}, 'UniformOutput', false))), 2); 
str_base_std = nanstd(full(cell2mat(cellfun(@(a) a(:, 1:19), {ss.unitTimeBStr}, 'UniformOutput', false))), 0, 2); 

for t = 1:length(ss)
    ss(t).utbStrZ = (full(ss(t).unitTimeBStr)-repmat(str_base_mean, 1, size(ss(t).unitTimeBStr, 2)))...
        ./repmat(str_base_std, 1, size(ss(t).unitTimeBStr, 2));
end

%% identify cells with significant silencing 
x = -975:50:3975; 
pre_stim = [];
post_stim = []; 
post_stim_min = []; 
breakthrough = []; 
reach = []; 
for t = 1:length(ss)
    if ~isempty(ss(t).utbStrStimAlignZ) && (isempty(ss(t).rStartRelToStim) || ss(t).rStartRelToStim>2000)
       pre_stim = [pre_stim, nanmean(ss(t).utbStrStimAlignZ(:, 1:20), 2)]; 
       post_stim_min = [post_stim, nanmin(ss(t).utbStrStimAlignZ(:, 21:40), [], 2)]; 
       post_stim = [post_stim, nanmean(ss(t).utbStrStimAlignZ(:, 21:40), 2)]; 
    end

    if ~isempty(ss(t).utbStrStimAlignZ) && ~isempty(ss(t).rStartRelToStim) && ss(t).rStartRelToStim<3500
        breakthrough = [breakthrough, nanmean(ss(t).utbStrStimAlignZ(:, x>ss(t).rStartRelToStim), 2)]; 
    end

    if strcmpi(ss(t).trialType, 'sp') && isempty(ss(t).tLaserStart)
        reach = [reach, nanmean(full(ss(t).utbStrZ), 21:40)]; 
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
%rez.gpfaCoeff = gpfaRezStr.faParams.L; 

% scatter_row_by_row_fillLogic([rez.magInh, rez.magBT, rez.magReach], [-1 3],  logical(ttestRez))
save(fullfile(filePath, 'magnitude_inhibition_breakThrough_str'), 'rez')






