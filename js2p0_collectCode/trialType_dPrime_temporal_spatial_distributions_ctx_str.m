%% get data
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez', 'dPrmTrjCell')
% Data comprise 9 sessions from 5 animals

isStrMat = cell2mat({dPrmRez(:).isStr})';
dPrmRezCtx = dPrmRez(~isStrMat);
dPrmRezStr = dPrmRez(isStrMat);

tType = cell2mat(cellfun(@(a) a(2), {dPrmRez.maxCoord_2dProj}, 'un', 0)); 
depth = [dPrmRez.depth]'; 
sigI = cell2mat({dPrmRez(:).sigI}');
sessionID = cellfun(@(a) a(1:11), {dPrmRez.cellId}, 'un', 0)'; 
mID = cellfun(@(a) a(1:4), {dPrmRez.cellId}, 'un', 0)'; 

%% Plot dPrm trajectories of direction, load, mixed encoding cells: STR
for tt = 1:8
    [mdPrm_str(tt,:),~,sdPrm_str(tt,:)] = meanstdsem(cell2mat(dPrmTrjCell([dPrmRez(:).isStr]',tt))); % # time bins: 101
end

min_sub_bool = true;
bins_to_plot = 25:75; 

dPrm_temporal.str_pure_load = min_max_norm(mdPrm_str([4, 8],:), min_sub_bool);  % pure load
dPrm_temporal.str_load_hi = min_max_norm(mdPrm_str(4,:), min_sub_bool);  % pure load - high 
dPrm_temporal.str_load_lo = min_max_norm(mdPrm_str(8,:), min_sub_bool);  % pure load - low 

dPrm_temporal.str_pure_dir = min_max_norm(mdPrm_str([2, 6],:), min_sub_bool);  % pure direction
dPrm_temporal.str_mixed = min_max_norm(mdPrm_str([1, 3, 5, 7],:), min_sub_bool);  % pure direction

% plot dir, load, mixed dPrm trajectory (distance from zero)
x = -0.5:.02:0.5;  
figure; hold on;
plot(x, dPrm_temporal.str_pure_dir(bins_to_plot))
plot(x, dPrm_temporal.str_pure_load(bins_to_plot))
plot(x, dPrm_temporal.str_mixed(bins_to_plot))
set(gca, 'TickDir', 'out')
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'dPrm_dist_from_zero_STR_dir_load_mixed'),'-dpdf','-vector')

% plot dir, load, mixed dPrm trajectory (distance from zero)
figure; hold on;
plot(dPrm_temporal.str_load_hi)
plot(dPrm_temporal.str_load_lo)
set(gca, 'TickDir', 'out')
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'dPrm_dist_from_zero_STR_hi_lo_load'),'-dpdf','-vector')

dPrm_temporal.str_temporal_rANOVA = rANOVA_dir_load_mixed(dPrmTrjCell([dPrmRez(:).isStr]', :));  

%% Plot dPrm spatial distribution of direction, load, mixed encoding cells: STR
min_sub_bool = false;
depth_str = depth(isStrMat); 
depth_str_bin = linspace(min(depth_str), max(depth_str), 19); %21); 
sessionID_str = sessionID(isStrMat); 
unique_sessionIDs_str = unique(sessionID_str); 
mID_str = mID(isStrMat); 
unique_mIDs_str = unique(mID_str); 

tt_str = tType(isStrMat); 
sig_str = sigI(isStrMat); 

depth_str_tt = cell(4, 8); 
for tt = 1:8
    % all counts
    depth_str_tt{1, tt} = depth_str((tt_str==tt)'); 
    [depth_str_tt{2, tt}, ~, depth_str_tt{3, tt}] = histcounts(depth_str_tt{1, tt}, depth_str_bin); % count
    depth_str_tt{4, tt} = full(smooth2a(depth_str_tt{2, tt}./sum(depth_str_tt{2, tt})*100, 0, 3)); % percent
    % significant counts
    depth_str_tt_sig{1, tt} = depth_str((tt_str==tt)' & sig_str);  
    [depth_str_tt_sig{2, tt}, ~, depth_str_tt_sig{3, tt}] = histcounts(depth_str_tt_sig{1, tt}, depth_str_bin); % count (significant cells only)
    depth_str_tt_sig{4, tt} = full(smooth2a(depth_str_tt_sig{2, tt}./sum(depth_str_tt_sig{2, tt})*100, 0, 3)); % percent (significant cells only)
    % per animal
    for mm = 1:length(unique(mID_str)) 
        temp_mI = cell2mat(cellfun(@(a) strcmpi(a, unique_mIDs_str{mm}), mID_str, 'un', 0)); 
        depth_str_tt_per_m{mm, tt} = depth_str((tt_str==tt)' & temp_mI); 
        depth_str_tt_per_m_sig{mm, tt} = depth_str((tt_str==tt)' & temp_mI & sig_str); 
    end
    % per session
    for ss = 1:length(unique(sessionID_str)) 
        temp_sI = cell2mat(cellfun(@(a) strcmpi(a, unique_sessionIDs_str{ss}), sessionID_str, 'un', 0)); 
        depth_str_tt_per_session{ss, tt} = depth_str((tt_str==tt)' & temp_sI); 
        depth_str_tt_per_session_sig{ss, tt} = depth_str((tt_str==tt)' & temp_sI & sig_str); 
    end
end

% 3 coarse bins STR counts
dPrm_spatial.str_pure_dir_counts_3bins  = sum(reshape(sum(cell2mat(depth_str_tt(2, [2, 6])')), [], 3)); % pure direction counts in 3 depth bins (superficial, intermediate, deep)
dPrm_spatial.str_pure_load_counts_3bins = sum(reshape(sum(cell2mat(depth_str_tt(2, [4, 8])')), [], 3)); % pure load counts in 3 depth bins
dPrm_spatial.str_mix_counts_3bins = sum(reshape(sum(cell2mat(depth_str_tt(2, [1, 3, 5, 7])')), [], 3)); % mixed counts in 3 depth bins
dPrm_spatial.str_pure_dir_load_counts_3bins = [dPrm_spatial.str_pure_dir_counts_3bins; dPrm_spatial.str_pure_load_counts_3bins]; 
myfisher(dPrm_spatial.str_pure_dir_load_counts_3bins); % counts without normalization 

% 3 coarse bins STR counts significant units only
dPrm_spatial.str_pure_dir_counts_3bins_sig  = sum(reshape(sum(cell2mat(depth_str_tt_sig(2, [2, 6])')), [], 3)); % pure direction counts in 3 depth bins (superficial, intermediate, deep)
dPrm_spatial.str_pure_load_counts_3bins_sig = sum(reshape(sum(cell2mat(depth_str_tt_sig(2, [4, 8])')), [], 3)); % pure load counts in 3 depth bins
dPrm_spatial.str_mix_counts_3bins_sig = sum(reshape(sum(cell2mat(depth_str_tt_sig(2, [1, 3, 5, 7])')), [], 3)); % mixed counts in 3 depth bins
dPrm_spatial.str_pure_dir_load_counts_3bins_sig = [dPrm_spatial.str_pure_dir_counts_3bins_sig; dPrm_spatial.str_pure_load_counts_3bins_sig]; 
myfisher(dPrm_spatial.str_pure_dir_load_counts_3bins_sig); % counts without normalization 

% 3 coarse bins STR fraction normalized to the MAX of the class (DO NOT USE THIS!)
dPrm_spatial.str_pure_dir = min_max_norm(cell2mat(depth_str_tt(4, [2, 6])'), min_sub_bool);  % pure direction
dPrm_spatial.str_pure_dir_fraction_3bins_norm = mean(reshape(dPrm_spatial.str_pure_dir, [], 3));  % pure direction
dPrm_spatial.str_pure_load = min_max_norm(cell2mat(depth_str_tt(4, [4, 8])'), min_sub_bool);   % pure load
dPrm_spatial.str_pure_load_fraction_3bins_norm = mean(reshape(dPrm_spatial.str_pure_load, [], 3));  % pure load
dPrm_spatial.str_mix = min_max_norm(cell2mat(depth_str_tt(4, [1, 3, 5, 7])'), min_sub_bool);  % mixed
dPrm_spatial.str_mix_fraction_3bins_norm = mean(reshape(dPrm_spatial.str_mix, [], 3));  % mixed
% 
% myfisher(round([dPrm_spatial.str_pure_dir_fraction_3bins_norm; dPrm_spatial.str_pure_load_fraction_3bins_norm])); % fraction without minmax normalization 

% 3 coarse bins STR fraction without normalization
dPrm_spatial.str_pure_dir_fraction_3bins = sum(reshape(mean(cell2mat(depth_str_tt(4, [2, 6])')), [], 3)); % pure direction fraction 3 bins  
dPrm_spatial.str_pure_load_fraction_3bins = sum(reshape(mean(cell2mat(depth_str_tt(4, [4, 8])')), [], 3)); % pure direction fraction 3 bins  
myfisher(round([dPrm_spatial.str_pure_dir_fraction_3bins; dPrm_spatial.str_pure_load_fraction_3bins])); % fraction without minmax normalization 

% 2 coarse bins STR fraction without normalization
dPrm_spatial.str_pure_dir_fraction_2bins = sum(reshape(mean(cell2mat(depth_str_tt(4, [2, 6])')), [], 2)); % pure direction fraction 2 bins  
dPrm_spatial.str_pure_load_fraction_2bins = sum(reshape(mean(cell2mat(depth_str_tt(4, [4, 8])')), [], 2)); % pure direction fraction 2 bins  
myfisher(round([dPrm_spatial.str_pure_dir_fraction_2bins; dPrm_spatial.str_pure_load_fraction_2bins])); % fraction without minmax normalization 

dPrm_spatial.str_pure_dir_fraction_2bins_sig = sum(reshape(mean(cell2mat(depth_str_tt_sig(4, [2, 6])')), [], 2)); % pure direction fraction 2 bins  
dPrm_spatial.str_pure_load_fraction_2bins_sig = sum(reshape(mean(cell2mat(depth_str_tt_sig(4, [4, 8])')), [], 2)); % pure direction fraction 2 bins  
myfisher(round([dPrm_spatial.str_pure_dir_fraction_2bins_sig; dPrm_spatial.str_pure_load_fraction_2bins_sig])); % fraction without minmax normalization 

% per animal mean depth paired t test 
dPrm_spatial.str_meanDepth_tt_per_m = cell2mat(cellfun(@mean, depth_str_tt_per_m, 'un', 0)); 
dPrm_spatial.str_meanDepth_tt_per_m_pure_dir = nanmean(dPrm_spatial.str_meanDepth_tt_per_m(:, [2, 6]), 2); 
dPrm_spatial.str_meanDepth_tt_per_m_pure_load = nanmean(dPrm_spatial.str_meanDepth_tt_per_m(:, [4, 8]), 2); 
[~, dPrm_spatial.str_meanDepth_per_m_paired_t_p, ~, dPrm_spatial.str_meanDepth_per_m_paired_t_stats] = ttest(dPrm_spatial.str_meanDepth_tt_per_m_pure_dir, ...
                                                                                                    dPrm_spatial.str_meanDepth_tt_per_m_pure_load); 
scatter_data(dPrm_spatial.str_meanDepth_tt_per_m_pure_dir, dPrm_spatial.str_meanDepth_tt_per_m_pure_load)
set(gca, 'YDir','reverse')

% per animal mean depth paired t test sig units only
dPrm_spatial.str_meanDepth_tt_per_m_sig = cell2mat(cellfun(@mean, depth_str_tt_per_m_sig, 'un', 0)); 
dPrm_spatial.str_meanDepth_tt_per_m_pure_dir_sig = nanmean(dPrm_spatial.str_meanDepth_tt_per_m_sig(:, [2, 6]), 2); 
dPrm_spatial.str_meanDepth_tt_per_m_pure_load_sig = nanmean(dPrm_spatial.str_meanDepth_tt_per_m_sig(:, [4, 8]), 2); 
[~, dPrm_spatial.str_meanDepth_per_m_sig_paired_t_p, ~, dPrm_spatial.str_meanDepth_per_m_sig_paired_t_stats] = ttest(dPrm_spatial.str_meanDepth_tt_per_m_pure_dir_sig, ...
                                                                                                    dPrm_spatial.str_meanDepth_tt_per_m_pure_load_sig); 
scatter_data(dPrm_spatial.str_meanDepth_tt_per_m_pure_dir_sig, dPrm_spatial.str_meanDepth_tt_per_m_pure_load_sig)
set(gca, 'YDir','reverse')

% per m delta (direction vs load) mean depth 
dPrm_spatial.str_meanDepth_dir_load_delta_depth = dPrm_spatial.str_meanDepth_tt_per_m_pure_dir - dPrm_spatial.str_meanDepth_tt_per_m_pure_load; 
dPrm_spatial.str_meanDepth_dir_load_delta_depth_sig = dPrm_spatial.str_meanDepth_tt_per_m_pure_dir_sig-dPrm_spatial.str_meanDepth_tt_per_m_pure_load_sig; 

% per session delta (direction vs load) mean depth
%dPrm_spatial.str_meanDepth_dir_load_delta_depth_per_session = dPrm_spatial.str_meanDepth_tt_per_session_pure_dir - dPrm_spatial.str_meanDepth_tt_per_session_pure_load; 
%dPrm_spatial.str_meanDepth_dir_load_delta_depth_sig_per_session = dPrm_spatial.str_meanDepth_tt_per_session_pure_dir_sig-dPrm_spatial.str_meanDepth_tt_per_session_pure_load_sig; 

% per session mean depth paired t test 
dPrm_spatial.str_meanDepth_tt_per_session = cell2mat(cellfun(@mean, depth_str_tt_per_session, 'un', 0)); 
dPrm_spatial.str_meanDepth_tt_per_session_pure_dir = nanmean(dPrm_spatial.str_meanDepth_tt_per_session(:, [2, 6]), 2); 
dPrm_spatial.str_meanDepth_tt_per_session_pure_load = nanmean(dPrm_spatial.str_meanDepth_tt_per_session(:, [4, 8]), 2); 
[~, dPrm_spatial.str_meanDepth_per_session_paired_t_p, ~, dPrm_spatial.str_meanDepth_per_session_paired_t_stats] = ttest(dPrm_spatial.str_meanDepth_tt_per_session_pure_dir, ...
                                                                                                    dPrm_spatial.str_meanDepth_tt_per_session_pure_load); 
scatter_data(dPrm_spatial.str_meanDepth_tt_per_session_pure_dir, dPrm_spatial.str_meanDepth_tt_per_session_pure_load)
set(gca, 'YDir','reverse')

% per session mean depth paired t test sig units only
dPrm_spatial.str_meanDepth_tt_per_session_sig = cell2mat(cellfun(@mean, depth_str_tt_per_session_sig, 'un', 0)); 
dPrm_spatial.str_meanDepth_tt_per_session_pure_dir_sig = nanmean(dPrm_spatial.str_meanDepth_tt_per_session_sig(:, [2, 6]), 2); 
dPrm_spatial.str_meanDepth_tt_per_session_pure_load_sig = nanmean(dPrm_spatial.str_meanDepth_tt_per_session_sig(:, [4, 8]), 2); 
[~, dPrm_spatial.str_meanDepth_per_session_sig_paired_t_p, ~, dPrm_spatial.str_meanDepth_per_session_sig_paired_t_stats] = ttest(dPrm_spatial.str_meanDepth_tt_per_session_pure_dir_sig, ...
                                                                                                    dPrm_spatial.str_meanDepth_tt_per_session_pure_load_sig); 
scatter_data(dPrm_spatial.str_meanDepth_tt_per_session_pure_dir_sig, dPrm_spatial.str_meanDepth_tt_per_session_pure_load_sig)
set(gca, 'YDir','reverse')

% rescale the str depths
min_str_depth = min([dPrmRez(isStrMat).depth]);
max_str_depth = max([dPrmRez(isStrMat).depth]); 
dPrm_spatial.rescaled_str_depth_dir = (cell2mat(depth_str_tt(1, [2, 6])')-min_str_depth)/max_str_depth; 
dPrm_spatial.rescaled_str_depth_load = (cell2mat(depth_str_tt(1, [4, 8])')-min_str_depth)/max_str_depth; 
dPrm_spatial.rescaled_str_depth_mix = (cell2mat(depth_str_tt(1, [1, 3, 5, 7])')-min_str_depth)/max_str_depth; 

[~, dPrm_spatial.rescaled_depth_ttest_str_dir_load_p, ~, dPrm_spatial.rescaled_depth_ttest_str_dir_load_stats] = ...
    ttest2(dPrm_spatial.rescaled_str_depth_dir, dPrm_spatial.rescaled_str_depth_load); 

% not rescaled
dPrm_spatial.str_depth_dir = cell2mat(depth_str_tt(1, [2, 6])'); 
dPrm_spatial.str_depth_load = cell2mat(depth_str_tt(1, [4, 8])'); 
dPrm_spatial.str_depth_mix = cell2mat(depth_str_tt(1, [1, 3, 5, 7])'); 
[~, dPrm_spatial.depth_ttest_str_dir_load_p, ~, dPrm_spatial.depth_ttest_str_dir_load_stats] = ...
    ttest2(cell2mat(depth_str_tt(1, [2, 6])'), cell2mat(depth_str_tt(1, [4, 8])')); 

% swarmchart of depth dir vs. load
figure; 
swarmchart([ones(1, length(dPrm_spatial.str_depth_dir)) 2*ones(1, length(dPrm_spatial.str_depth_load))], ...
    [dPrm_spatial.str_depth_dir', dPrm_spatial.str_depth_load'])


figure; hold on; 
plot(depth_str_bin(1:end-1), dPrm_spatial.str_pure_dir, '-O')
plot(depth_str_bin(1:end-1), dPrm_spatial.str_pure_load, '-O')
plot(depth_str_bin(1:end-1), dPrm_spatial.str_mix, '-O')
set(gca, 'TickDir', 'out')
ylim([0.2, 1])
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'dPrm_spatial_distribute_STR'),'-dpdf','-vector')

% corr spatial distribution of direction- and load- tuned neurons
dPrm_spatial.depth_str_dir = mean(cell2mat(depth_str_tt(4, [2, 6])'))'; 
dPrm_spatial.depth_str_load = mean(cell2mat(depth_str_tt(4, [4, 8])'))'; 
dPrm_spatial.depth_str_mix = mean(cell2mat(depth_str_tt(4, [1, 3, 5, 7])'))'; 

[dPrm_spatial.corr_Rho_str_dir_load, dPrm_spatial.corr_P_str_dir_load] = corr(dPrm_spatial.depth_str_dir, dPrm_spatial.depth_str_load); 
[dPrm_spatial.corr_Rho_str_mix_load, dPrm_spatial.corr_P_str_mix_load] = corr(dPrm_spatial.depth_str_mix, dPrm_spatial.depth_str_load); 
[dPrm_spatial.corr_Rho_str_dir_mix, dPrm_spatial.corr_P_str_dir_mix] = corr(dPrm_spatial.depth_str_dir, dPrm_spatial.depth_str_mix); 

%% Plot dPrm trajectories of direction, load, mixed encoding cells: CTX
for tt = 1:8
    [mdPrm_ctx(tt,:),~,sdPrm_ctx(tt,:)] = meanstdsem(cell2mat(dPrmTrjCell(~[dPrmRez(:).isStr]',tt)));
end

min_sub_bool = true; 
dPrm_temporal.ctx_pure_load = min_max_norm(mdPrm_ctx([4, 8],:), min_sub_bool);  % pure load
dPrm_temporal.ctx_load_hi = min_max_norm(mdPrm_ctx(4,:), min_sub_bool);  % pure load - high 
dPrm_temporal.ctx_load_lo = min_max_norm(mdPrm_ctx(8,:), min_sub_bool);  % pure load - low 

dPrm_temporal.ctx_pure_dir = min_max_norm(mdPrm_ctx([2, 6],:), min_sub_bool);  % pure direction
dPrm_temporal.ctx_mixed = min_max_norm(mdPrm_ctx([1, 3, 5, 7],:), min_sub_bool);  % pure direction

% plot dir, load, mixed dPrm trajectory (distance from zero)
figure; hold on;
plot(x, dPrm_temporal.ctx_pure_dir(bins_to_plot))
plot(x, dPrm_temporal.ctx_pure_load(bins_to_plot))
plot(x, dPrm_temporal.ctx_mixed(bins_to_plot))
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'dPrm_dist_from_zero_CTX_dir_load_mixed'),'-dpdf','-vector')

% plot dir, load, mixed dPrm trajectory (distance from zero)
figure; hold on;
plot(dPrm_temporal.ctx_load_hi)
plot(dPrm_temporal.ctx_load_lo)
set(gca, 'TickDir', 'out')
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'dPrm_dist_from_zero_CTX_hi_lo_load'),'-dpdf','-vector')

dPrm_temporal.ctx_temporal_rANOVA = rANOVA_dir_load_mixed(dPrmTrjCell(~[dPrmRez(:).isStr]', :)); 

%% Plot dPrm spatial dictxibution of direction, load, mixed encoding cells: CTX
min_sub_bool = false;
depth_ctx = depth(~isStrMat); 
depth_ctx_bin = linspace(min(depth_ctx), max(depth_ctx), 19); %21); 
sessionID_ctx = sessionID(~isStrMat); 
unique_sessionIDs_ctx = unique(sessionID_ctx); 
mID_ctx = mID(~isStrMat); 
unique_mIDs_ctx = unique(mID_ctx); 

tt_ctx = tType(~isStrMat); 
sig_ctx = sigI(~isStrMat); 

depth_ctx_tt = cell(3, 8); 
for tt = 1:8
    % all counts
    depth_ctx_tt{1, tt} = depth_ctx((tt_ctx==tt)'); 
    [depth_ctx_tt{2, tt}, ~, depth_ctx_tt{3, tt}] = histcounts(depth_ctx_tt{1, tt}, depth_ctx_bin); % count
    depth_ctx_tt{4, tt} = full(smooth2a(depth_ctx_tt{2, tt}./sum(depth_ctx_tt{2, tt})*100, 0, 3)); % percent
    % significant counts
    depth_ctx_tt_sig{1, tt} = depth_ctx((tt_ctx==tt)' & sig_ctx);  
    [depth_ctx_tt_sig{2, tt}, ~, depth_ctx_tt_sig{3, tt}] = histcounts(depth_ctx_tt_sig{1, tt}, depth_ctx_bin); % count (significant cells only)
    depth_ctx_tt_sig{4, tt} = full(smooth2a(depth_ctx_tt_sig{2, tt}./sum(depth_ctx_tt_sig{2, tt})*100, 0, 3)); % percent (significant cells only)
    % per animal
    for mm = 1:length(unique(mID_ctx)) 
        temp_mI = cell2mat(cellfun(@(a) strcmpi(a, unique_mIDs_ctx{mm}), mID_ctx, 'un', 0)); 
        depth_ctx_tt_per_m{mm, tt} = depth_ctx((tt_ctx==tt)' & temp_mI); 
        depth_ctx_tt_per_m_sig{mm, tt} = depth_ctx((tt_ctx==tt)' & temp_mI & sig_ctx); 
    end
    % per session
    for ss = 1:length(unique(sessionID_ctx)) 
        temp_sI = cell2mat(cellfun(@(a) strcmpi(a, unique_sessionIDs_ctx{ss}), sessionID_ctx, 'un', 0)); 
        depth_ctx_tt_per_session{ss, tt} = depth_ctx((tt_ctx==tt)' & temp_sI); 
        depth_ctx_tt_per_session_sig{ss, tt} = depth_ctx((tt_ctx==tt)' & temp_sI & sig_ctx); 
    end
end

% 3 coarse bins CTX counts
dPrm_spatial.ctx_pure_dir_counts_3bins  = sum(reshape(sum(cell2mat(depth_ctx_tt(2, [2, 6])')), [], 3)); % pure direction counts in 3 depth bins (superficial, intermediate, deep)
dPrm_spatial.ctx_pure_load_counts_3bins = sum(reshape(sum(cell2mat(depth_ctx_tt(2, [4, 8])')), [], 3)); % pure load counts in 3 depth bins
dPrm_spatial.ctx_mix_counts_3bins = sum(reshape(sum(cell2mat(depth_ctx_tt(2, [1, 3, 5, 7])')), [], 3)); % mixed counts in 3 depth bins
myfisher([dPrm_spatial.ctx_pure_dir_counts_3bins; dPrm_spatial.ctx_pure_load_counts_3bins]); % counts without normalization 

% 3 coarse bins CTX counts significant units only
dPrm_spatial.ctx_pure_dir_counts_3bins_sig  = sum(reshape(sum(cell2mat(depth_ctx_tt_sig(2, [2, 6])')), [], 3)); % pure direction counts in 3 depth bins (superficial, intermediate, deep)
dPrm_spatial.ctx_pure_load_counts_3bins_sig = sum(reshape(sum(cell2mat(depth_ctx_tt_sig(2, [4, 8])')), [], 3)); % pure load counts in 3 depth bins
dPrm_spatial.ctx_mix_counts_3bins_sig = sum(reshape(sum(cell2mat(depth_ctx_tt_sig(2, [1, 3, 5, 7])')), [], 3)); % mixed counts in 3 depth bins
dPrm_spatial.ctx_pure_dir_load_counts_3bins_sig = [dPrm_spatial.ctx_pure_dir_counts_3bins_sig; dPrm_spatial.ctx_pure_load_counts_3bins_sig]; 
myfisher(dPrm_spatial.ctx_pure_dir_load_counts_3bins_sig); % counts without normalization 

% 3 coarse bins CTX fraction normalized to the MAX of the class (DO NOT USE THIS!)
dPrm_spatial.ctx_pure_dir = min_max_norm(cell2mat(depth_ctx_tt(4, [2, 6])'), min_sub_bool);  % pure direction
dPrm_spatial.ctx_pure_dir_fraction_3bins_norm = mean(reshape(dPrm_spatial.ctx_pure_dir, [], 3));  % pure direction
dPrm_spatial.ctx_pure_load = min_max_norm(cell2mat(depth_ctx_tt(4, [4, 8])'), min_sub_bool);   % pure load
dPrm_spatial.ctx_pure_load_fraction_3bins_norm = mean(reshape(dPrm_spatial.ctx_pure_load, [], 3));  % pure load
dPrm_spatial.ctx_mix = min_max_norm(cell2mat(depth_ctx_tt(4, [1, 3, 5, 7])'), min_sub_bool);  % mixed
dPrm_spatial.ctx_mix_fraction_3bins_norm = mean(reshape(dPrm_spatial.ctx_mix, [], 3));  % mixed
%
% myfisher(round([dPrm_spatial.ctx_pure_dir_fraction_3bins_norm; dPrm_spatial.ctx_pure_load_fraction_3bins_norm])); % fraction without minmax normalization 

% 3 coarse bins CTX fraction without normalization
dPrm_spatial.ctx_pure_dir_fraction_3bins = sum(reshape(mean(cell2mat(depth_ctx_tt(4, [2, 6])')), [], 3)); % pure direction fraction 3 bins  
dPrm_spatial.ctx_pure_load_fraction_3bins = sum(reshape(mean(cell2mat(depth_ctx_tt(4, [4, 8])')), [], 3)); % pure direction fraction 3 bins  
myfisher(round([dPrm_spatial.ctx_pure_dir_fraction_3bins; dPrm_spatial.ctx_pure_load_fraction_3bins])); % fraction without minmax normalization 

% 2 coarse bins CTX fraction without normalization
dPrm_spatial.ctx_pure_dir_fraction_2bins = sum(reshape(mean(cell2mat(depth_ctx_tt(4, [2, 6])')), [], 2)); % pure direction fraction 2 bins  
dPrm_spatial.ctx_pure_load_fraction_2bins = sum(reshape(mean(cell2mat(depth_ctx_tt(4, [4, 8])')), [], 2)); % pure direction fraction 2 bins  
myfisher(round([dPrm_spatial.ctx_pure_dir_fraction_2bins; dPrm_spatial.ctx_pure_load_fraction_2bins])); % fraction without minmax normalization 

dPrm_spatial.ctx_pure_dir_fraction_2bins_sig = sum(reshape(mean(cell2mat(depth_ctx_tt_sig(4, [2, 6])')), [], 2)); % pure direction fraction 2 bins  
dPrm_spatial.ctx_pure_load_fraction_2bins_sig = sum(reshape(mean(cell2mat(depth_ctx_tt_sig(4, [4, 8])')), [], 2)); % pure direction fraction 2 bins  
myfisher(round([dPrm_spatial.ctx_pure_dir_fraction_2bins_sig; dPrm_spatial.ctx_pure_load_fraction_2bins_sig])); % fraction without minmax normalization 

% per animal mean depth paired t test 
dPrm_spatial.ctx_meanDepth_tt_per_m = cell2mat(cellfun(@mean, depth_ctx_tt_per_m, 'un', 0)); 
dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_m(:, [2, 6]), 2); 
dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_m(:, [4, 8]), 2); 
[~, dPrm_spatial.ctx_meanDepth_per_m_paired_t_p, ~, dPrm_spatial.ctx_meanDepth_per_m_paired_t_stats] = ttest(dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir, ...
                                                                                                    dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load); 
scatter_data(dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir, dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load)
set(gca, 'YDir','reverse')

% per animal mean depth paired t test sig units only
dPrm_spatial.ctx_meanDepth_tt_per_m_sig = cell2mat(cellfun(@mean, depth_ctx_tt_per_m_sig, 'un', 0)); 
dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir_sig = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_m_sig(:, [2, 6]), 2); 
dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load_sig = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_m_sig(:, [4, 8]), 2); 
[~, dPrm_spatial.ctx_meanDepth_per_m_sig_paired_t_p, ~, dPrm_spatial.ctx_meanDepth_per_m_sig_paired_t_stats] = ttest(dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir_sig, ...
                                                                                                    dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load_sig); 
scatter_data(dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir_sig, dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load_sig)
set(gca, 'YDir','reverse')

% per m delta (direction vs load) mean depth 
dPrm_spatial.ctx_meanDepth_dir_load_delta_depth = dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir - dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load; 
dPrm_spatial.ctx_meanDepth_dir_load_delta_depth_sig = dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir_sig-dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load_sig; 

% per session mean depth paired t test 
dPrm_spatial.ctx_meanDepth_tt_per_session = cell2mat(cellfun(@mean, depth_ctx_tt_per_session, 'un', 0)); 
dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_session(:, [2, 6]), 2); 
dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_session(:, [4, 8]), 2); 
[~, dPrm_spatial.ctx_meanDepth_per_session_paired_t_p, ~, dPrm_spatial.ctx_meanDepth_per_session_paired_t_stats] = ttest(dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir, ...
                                                                                                    dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load); 
scatter_data(dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir, dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load)
set(gca, 'YDir','reverse')

% per session mean depth paired t test sig units only
dPrm_spatial.ctx_meanDepth_tt_per_session_sig = cell2mat(cellfun(@mean, depth_ctx_tt_per_session_sig, 'un', 0)); 
dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir_sig = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_session_sig(:, [2, 6]), 2); 
dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load_sig = nanmean(dPrm_spatial.ctx_meanDepth_tt_per_session_sig(:, [4, 8]), 2); 
[~, dPrm_spatial.ctx_meanDepth_per_session_sig_paired_t_p, ~, dPrm_spatial.ctx_meanDepth_per_session_sig_paired_t_stats] = ttest(dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir_sig, ...
                                                                                                    dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load_sig); 
scatter_data(dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir_sig, dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load_sig)
set(gca, 'YDir','reverse')

% per m delta (direction vs load) mean depth 
dPrm_spatial.ctx_meanDepth_dir_load_delta_depth = dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir - dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load; 
dPrm_spatial.ctx_meanDepth_dir_load_delta_depth_sig = dPrm_spatial.ctx_meanDepth_tt_per_m_pure_dir_sig-dPrm_spatial.ctx_meanDepth_tt_per_m_pure_load_sig; 

% per session delta (direction vs load) mean depth
dPrm_spatial.ctx_meanDepth_dir_load_delta_depth_per_session = dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir - dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load; 
dPrm_spatial.ctx_meanDepth_dir_load_delta_depth_sig_per_session = dPrm_spatial.ctx_meanDepth_tt_per_session_pure_dir_sig-dPrm_spatial.ctx_meanDepth_tt_per_session_pure_load_sig; 

% rescale the ctx depths
min_ctx_depth = min([dPrmRez(~isStrMat).depth]);
max_ctx_depth = max([dPrmRez(~isStrMat).depth]); 
dPrm_spatial.rescaled_ctx_depth_dir = (cell2mat(depth_ctx_tt(1, [2, 6])')-min_ctx_depth)/max_ctx_depth; 
dPrm_spatial.rescaled_ctx_depth_load = (cell2mat(depth_ctx_tt(1, [4, 8])')-min_ctx_depth)/max_ctx_depth; 
dPrm_spatial.rescaled_ctx_depth_mix = (cell2mat(depth_ctx_tt(1, [1, 3, 5, 7])')-min_ctx_depth)/max_ctx_depth; 

% rescaled
[~, dPrm_spatial.rescaled_depth_ttest_ctx_dir_load_p, ~, dPrm_spatial.rescaled_depth_ttest_ctx_dir_load_stats] = ...
    ttest2(dPrm_spatial.rescaled_ctx_depth_dir, dPrm_spatial.rescaled_ctx_depth_load); 

% not rescaled
dPrm_spatial.ctx_depth_dir = cell2mat(depth_ctx_tt(1, [2, 6])'); 
dPrm_spatial.ctx_depth_load = cell2mat(depth_ctx_tt(1, [4, 8])'); 
dPrm_spatial.ctx_depth_mix = cell2mat(depth_ctx_tt(1, [1, 3, 5, 7])'); 
[~, dPrm_spatial.depth_ttest_ctx_dir_load_p, ~, dPrm_spatial.depth_ttest_ctx_dir_load_stats] = ...
    ttest2(cell2mat(depth_ctx_tt(1, [2, 6])'), cell2mat(depth_ctx_tt(1, [4, 8])')); 

% swarmchart of depth dir vs. load
figure; 
swarmchart([ones(1, length(dPrm_spatial.ctx_depth_dir)) 2*ones(1, length(dPrm_spatial.ctx_depth_load))], ...
    [dPrm_spatial.ctx_depth_dir', dPrm_spatial.ctx_depth_load'])

% more plots
figure; hold on; 
plot(depth_ctx_bin(1:end-1), dPrm_spatial.ctx_pure_dir, '-O')
plot(depth_ctx_bin(1:end-1), dPrm_spatial.ctx_pure_load, '-O')
plot(depth_ctx_bin(1:end-1), dPrm_spatial.ctx_mix, '-O')
set(gca, 'TickDir', 'out')
ylim([0.2, 1])
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'dPrm_spatial_distribute_CTX'),'-dpdf','-vector')

% corr spatial distribution of direction- and load- tuned neurons
dPrm_spatial.depth_ctx_dir = mean(cell2mat(depth_ctx_tt(4, [2, 6])'))'; 
dPrm_spatial.depth_ctx_load = mean(cell2mat(depth_ctx_tt(4, [4, 8])'))'; 
dPrm_spatial.depth_ctx_mix = mean(cell2mat(depth_ctx_tt(4, [1, 3, 5, 7])'))'; 

[dPrm_spatial.corr_Rho_ctx_dir_load, dPrm_spatial.corr_P_ctx_dir_load] = corr(dPrm_spatial.depth_ctx_dir, dPrm_spatial.depth_ctx_load); 
[dPrm_spatial.corr_Rho_ctx_mix_load, dPrm_spatial.corr_P_ctx_mix_load] = corr(dPrm_spatial.depth_ctx_mix, dPrm_spatial.depth_ctx_load); 
[dPrm_spatial.corr_Rho_ctx_dir_mix, dPrm_spatial.corr_P_ctx_dir_mix] = corr(dPrm_spatial.depth_ctx_dir, dPrm_spatial.depth_ctx_mix); 

%% SAVE dPrm_spatial and dPrm_temporal for further analysis and plotting
save(fullfile('/Volumes/Extreme SSD/js2p0/collectData','dPrm_temporal_spatial_rez_CtxStr'), 'dPrm_temporal', 'dPrm_spatial')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x_n = min_max_norm(x, min_sub)
x_m = nanmean(x, 1); 

if min_sub == true 
    x_m_min = min(x_m); 
    x_m = x_m - x_m_min; 
end

x_m_max = max(x_m);
x_n = x_m/x_m_max; 
end

function rANOVA_rez = rANOVA_dir_load_mixed(dPrmTrjC) 
% dPrmTrjC = dPrmTrjCell([dPrmRez(:).isStr]', :); 
% encoding types 
dir_tt = [2, 6]; 
load_tt = [4, 8]; 
mix_tt = [1, 3, 5, 7]; 

% time bins to be used
bins = 21:80; 

mm_dPrmTrjC = cell(1, size(dPrmTrjC, 2)); 

for t = 1:size(dPrmTrjC, 2)
    valI = ~cell2mat(cellfun(@isempty, dPrmTrjC(:, t), 'un', 0)); 
    
    mm_dPrmTrjC{1, t} = cell2mat(cellfun(@(a) min_max_norm(a, true), dPrmTrjC(valI, t), 'un', 0)); 
end

%avg_mm_dPrmTrjC = cellfun(@mean, mm_dPrmTrjC, 'un', 0); 

mm_dPrmTrjC_dir = cell2mat(cellfun(@(a) a(:, bins), mm_dPrmTrjC(dir_tt)', 'un', 0)); 
mm_dPrmTrjC_load = cell2mat(cellfun(@(a) a(:, bins), mm_dPrmTrjC(load_tt)', 'un', 0)); 
mm_dPrmTrjC_mix = cell2mat(cellfun(@(a) a(:, bins), mm_dPrmTrjC(mix_tt)', 'un', 0)); 

anova_mat = [mm_dPrmTrjC_dir; mm_dPrmTrjC_load; mm_dPrmTrjC_mix]; 

% group variable
group_var_dir = cell(size(mm_dPrmTrjC_dir, 1), 1); 
[group_var_dir{:}] = deal('dir'); 

group_var_load = cell(size(mm_dPrmTrjC_load, 1), 1); 
[group_var_load{:}] = deal('load'); 

group_var_mix = cell(size(mm_dPrmTrjC_mix, 1), 1); 
[group_var_mix{:}] = deal('mix');

group_var = [group_var_dir; group_var_load; group_var_mix]; 

% variable names 
var_names = cell(1, 1+size(anova_mat, 2)); 
for jj = 1:length(var_names)
    if jj == 1
        var_names{jj} = 'G'; 
    else
        var_names{jj} = strcat('B', num2str(jj-1)); 
    end
end

% build table
rm_tab = cell2table([group_var, num2cell(anova_mat)], 'VariableNames', var_names); 

Block = table((1:size(anova_mat,2))','VariableNames',{'Block'}); % Block
% Fit repetitive model to data
rm_design = ['B1-', strcat('B', num2str(size(anova_mat,2))), ' ', '~', ' ', 'G'];  % rm design string

rANOVA_rez.rm = fitrm(rm_tab, rm_design, 'WithinDesign', Block); % repeated measures model fit for RT data
rANOVA_rez.ranova_rez = ranova(rm); % repeated measures ANOVA on RT (within- and within*group interaction)
rANOVA_rez.multiCompBlock = multcompare(rm, 'G', 'By', 'Block');

end

function scatter_data(dataset1, varargin)
numb_sets = length(varargin) + 1; 

figure; hold on; 

n = numel(dataset1); 

for k = 1:numb_sets
    if k == 1
        x_jitter = (rand(length(dataset1), 1)-0.5)/10; 
        x = k*ones(length(dataset1), 1) + x_jitter; 
        y = dataset1; 
        scatter(x, y, 'k'); 
        pre_data = {x, y}; 
    else
        assert(n==numel(varargin{k-1}))
        x_jitter = (rand(length(dataset1), 1)-0.5)/10; 
        x = k*ones(length(varargin{k-1}), 1) + x_jitter;
        y = varargin{k-1}; 
        scatter(x, y, 'k'); 
        % draw line plot
        for kk = 1:n
            plot([pre_data{1}(kk), x(kk)], [pre_data{2}(kk), y(kk)], 'k:')
        end
        pre_data = {x, y}; 
    end
end
hold off; 
xlim([0.5 n+0.5])

end 
