
filePath_reach = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/rezKFdecodeHTrjCtxStrPosVel_reach_new_WR40_082019.mat'; 
filePath_pull = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/rezKFdecodeHTrjCtxStrPosVel_pull_new_WR40_082019.mat'; 

figure_name_header = 'WR40_082019'; 

%% load REACH data
load(fullfile(filePath_reach), 's', 'rez_reach')
reachState = s.dat.stateR(:, 1); 

% median subtraction 
valStateI = ~cell2mat(cellfun(@isempty, s.dat.stateR, 'un', 0)); 
valStates = s.dat.stateR(valStateI); 
valStates_init_med = median(cell2mat(cellfun(@(a) a(1:3, 1), valStates, 'un', 0)'), 2); 

reachState = reachState(valStateI(:, 1));
reachState = cellfun(@(a) a(1:3, :), reachState, 'un', 0); 
reachState_medSub = cellfun(@(a) a-valStates_init_med, reachState, 'un', 0); 

reachState_ctx = rez_reach{1, 7}.rpr_est_ctx(valStateI(:, 1), 1); % representative trajectories estimated on MCtx activity
reachState_str = rez_reach{1, 7}.rpr_est_str(valStateI(:, 1), 1); % representative trajectories estimated on Str activity
reachState_cg = rez_reach{1, 7}.rpr_est_cg(valStateI(:, 1), 1); % representative trajectories estimated on Cg activity
clearvars s 

%% load PULL data
load(fullfile(filePath_pull), 's', 'rez_pull')
pullState = s.dat.stateP(:, 1); 

% median subtraction 
valStateI = ~cell2mat(cellfun(@isempty, s.dat.stateP, 'un', 0)); 
valStates_pull = s.dat.stateP(valStateI); 
valStates_pull_init_med = median(cell2mat(cellfun(@(a) a(1:3, 1), valStates_pull, 'un', 0)'), 2); 

pullState = pullState(valStateI(:, 1));
pullState = cellfun(@(a) a(1:3, :), pullState, 'un', 0); 
pullState_medSub = cellfun(@(a) a-valStates_init_med, pullState, 'un', 0); 
reachPull_init_diff = valStates_init_med - valStates_pull_init_med; 

pullState_ctx_wo_correction = rez_pull{1, 7}.rpr_est_ctx(valStateI(:, 1), 1); % representative trajectories estimated on MCtx activity
pullState_str_wo_correction = rez_pull{1, 7}.rpr_est_str(valStateI(:, 1), 1); % representative trajectories estimated on Str activity
pullState_cg_wo_correction = rez_pull{1, 7}.rpr_est_cg(valStateI(:, 1), 1); % representative trajectories estimated on Cg activity
clearvars s 

pullState_ctx = cellfun(@(a) a-reachPull_init_diff, pullState_ctx_wo_correction, 'UniformOutput', false); 
pullState_str = cellfun(@(a) a-reachPull_init_diff, pullState_str_wo_correction, 'UniformOutput', false); 
pullState_cg = cellfun(@(a) a-reachPull_init_diff, pullState_cg_wo_correction, 'UniformOutput', false); 

%% combine reach and pull
reachPullState = cellfun(@(a, b) [a, b], reachState_medSub, pullState_medSub, 'UniformOutput', false); 
reachPullState_ctx = cellfun(@(a, b) [a, b], reachState_ctx, pullState_ctx, 'UniformOutput', false); 
reachPullState_str = cellfun(@(a, b) [a, b], reachState_str, pullState_str, 'UniformOutput', false); 
reachPullState_cg = cellfun(@(a, b) [a, b], reachState_cg, pullState_cg, 'UniformOutput', false); 

%% correlation for trial sorting
corr_reachPull = cellfun(@(a, b) corrcoef(a(2, :), b(2, :)), reachPullState, reachPullState_ctx, 'un', 0); 
corr_reachPull = cell2mat(cellfun(@(a) a(1, 2), corr_reachPull, 'un', 0)); 
corr_reachPull_sort = sortrows([corr_reachPull, (1:length(corr_reachPull))', cell2mat(cellfun(@length, reachState_medSub, 'un', 0))], -1); 

%% Y-axis (front - back)
figure; hold on; 
x_init = 5; 
for t = [1:2, 4:11]
    trialI = corr_reachPull_sort(t, 2); 
    reachEnd = corr_reachPull_sort(t, 3); 

    % get state trj
    rpState = smooth2a(reachPullState{trialI}(2, :), 0, 1); 
    rpState_ctx = smooth2a(reachPullState_ctx{trialI}(2, :), 0, 1); 
    rpState_str = smooth2a(reachPullState_str{trialI}(2, :), 0, 1); 
    rpState_cg = smooth2a(reachPullState_cg{trialI}(2, :), 0, 1); 
    
    % get x coordinates
    x = x_init:x_init+length(rpState_cg)-1; 
    
    % patch pull phase
    patch([x_init+reachEnd, x_init+length(rpState), x_init+length(rpState), x_init+reachEnd], ...
              [-4 -4 14 14], [235, 176, 32]./255, 'LineStyle','none')
    alpha(0.2)

    alpha(0.7)
    % plot trajectories
    plot(x, -rpState, "k"); 
    plot(x, -rpState_ctx, "Color", [27 117 187]./255); 
    plot(x, -rpState_str, "Color", [0 166 156]./255); 
    plot(x, -rpState_cg, "Color", [158 31 99]./255); 
    
    scatter(x(1), -rpState(1), 150, "k", "filled")
    scatter(x(1), -rpState_ctx(1), 150, [27 117 187]./255, "filled")
    scatter(x(1), -rpState_str(1), 150, [0 166 156]./255, "filled")
    scatter(x(1), -rpState_cg(1), 150, [158 31 99]./255, "filled")

    x_init = x(end) + 5; 
end
xlim([0 300])
ylim([-4 14])
set(gca, 'TickDir', 'out')

print(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'KF_representative trajectories_Y_axis'), '-dpdf', '-vector')


%% X-axis (front - back)
figure; hold on; 
x_init = 5; 
for t = [1:2, 4:11]
    trialI = corr_reachPull_sort(t, 2); 
    reachEnd = corr_reachPull_sort(t, 3); 

    % get state trj
    rpState = smooth2a(reachPullState{trialI}(1, :), 0, 1); 
    rpState_ctx = smooth2a(reachPullState_ctx{trialI}(1, :), 0, 1); 
    rpState_str = smooth2a(reachPullState_str{trialI}(1, :), 0, 1); 
    rpState_cg = smooth2a(reachPullState_cg{trialI}(1, :), 0, 1); 
    
    % get x coordinates
    x = x_init:x_init+length(rpState_cg)-1; 
    
    % patch pull phase
    patch([x_init+reachEnd, x_init+length(rpState), x_init+length(rpState), x_init+reachEnd], ...
              [4 4 -10 -10], [235, 176, 32]./255, 'LineStyle','none')
    alpha(0.2)

    alpha(0.7)
    % plot trajectories
    plot(x, -rpState, "k"); 
    plot(x, -rpState_ctx, "Color", [27 117 187]./255); 
    plot(x, -rpState_str, "Color", [0 166 156]./255); 
    plot(x, -rpState_cg, "Color", [158 31 99]./255); 
    
    scatter(x(1), -rpState(1), 150, "k", "filled")
    scatter(x(1), -rpState_ctx(1), 150, [27 117 187]./255, "filled")
    scatter(x(1), -rpState_str(1), 150, [0 166 156]./255, "filled")
    scatter(x(1), -rpState_cg(1), 150, [158 31 99]./255, "filled")

    x_init = x(end) + 5; 
end
xlim([0 300])
ylim([-10 3])
set(gca, 'TickDir', 'out')

print(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'KF_representative trajectories_X_axis'), '-dpdf', '-vector')

%% Z-axis (front - back)
figure; hold on; 
x_init = 5; 
for t = [1:2, 4:11]
    trialI = corr_reachPull_sort(t, 2); 
    reachEnd = corr_reachPull_sort(t, 3); 

    % get state trj
    rpState = smooth2a(reachPullState{trialI}(3, :), 0, 1); 
    rpState_ctx = smooth2a(reachPullState_ctx{trialI}(3, :), 0, 1); 
    rpState_str = smooth2a(reachPullState_str{trialI}(3, :), 0, 1); 
    rpState_cg = smooth2a(reachPullState_cg{trialI}(3, :), 0, 1); 
    
    % get x coordinates
    x = x_init:x_init+length(rpState_cg)-1; 
    
    % patch pull phase
    patch([x_init+reachEnd, x_init+length(rpState), x_init+length(rpState), x_init+reachEnd], ...
              [-4 -4 9 9], [235, 176, 32]./255, 'LineStyle','none')
    alpha(0.2)

    alpha(0.7)
    % plot trajectories
    plot(x, rpState, "k"); 
    plot(x, rpState_ctx, "Color", [27 117 187]./255); 
    plot(x, rpState_str, "Color", [0 166 156]./255); 
    plot(x, rpState_cg, "Color", [158 31 99]./255); 
    
    scatter(x(1), rpState(1), 150, "k", "filled")
    scatter(x(1), rpState_ctx(1), 150, [27 117 187]./255, "filled")
    scatter(x(1), rpState_str(1), 150, [0 166 156]./255, "filled")
    scatter(x(1), rpState_cg(1), 150, [158 31 99]./255, "filled")

    x_init = x(end) + 5; 
end
xlim([0 300])
ylim([-4 9])
set(gca, 'TickDir', 'out')

print(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'KF_representative trajectories_Z_axis'), '-dpdf', '-vector')



