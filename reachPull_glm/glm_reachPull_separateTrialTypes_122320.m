% Fit GLM a version 

clc; clearvars; close all;

%% 1. Load data
filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'; 
cd(filePath)

% neural and behavioral data
spkDir = dir('binSpkCountSTRCTX*');
load(fullfile(spkDir(1).folder, spkDir(1).name),'spkTimesCell','jkvt','meta','rStartToPull')
S=rStartToPull; 

% task parameter
DT = 10; % 10 ms (width of timeBin)
N_MOVE = 10; % reach start, n bumps
RANGE_MOVE = [0 3000]; 
BIAS_MOVE = 5;
N_TASK = 10; % reward
RANGE_TASK = [0 3000]; % ms (0 to 3000 ms relative to reward delivery) 
BIAS_TASK = 8;
N_HANDVEL = 10; 

% code parameter
PLOT = true;
CALC_PRM = true;

%% 2. Behavioral data 
% detect events and continuous joystick pull trajectory  
rwdTimes = cell2mat(cellfun(@(a) a(~isnan(a)), {jkvt(:).rewardT}, 'un', 0)); 
%sessionDur = round(str2double(meta.fileTimeSecs)*1000); % in ms 
sessionDur = rwdTimes(end)+5*1000; % cut 5 seconds after the last reward
session_range = [1 sessionDur]; 
clip = @(t) t(session_range(1) <= t & t <= session_range(2)) - session_range(1);
time_bin = (0:DT:sessionDur)';
n_bin = length(time_bin) - 1;

% trial indices
[posTqC,posTqTypes] = posTrqCombinations( jkvt ); 
rStartI = cellfun(@(c) ~isempty(c), {jkvt(:).rStartToPull});
pStartI = cellfun(@(c) ~isempty(c), {jkvt(:).pullStarts}); 
leI = cell2mat(cellfun(@(a) contains(a,'p1'), posTqC, 'un', 0 )) &pStartI'; 
riI = cell2mat(cellfun(@(a) contains(a,'p2'), posTqC, 'un', 0 )) &pStartI';  
loI = cell2mat(cellfun(@(a) contains(a,'t1'), posTqC, 'un', 0 )) &pStartI'; 
hiI = cell2mat(cellfun(@(a) contains(a,'t2'), posTqC, 'un', 0 )) &pStartI';  
rwI   = [jkvt(:).rewarded]'; 

% if reach Start is missing, put estimated values based on the pullStarts with (-200 ms offset)
for tt = find(~rStartI&pStartI)
    jkvt(tt).rStartToPull = jkvt(tt).pullStarts-200; 
end
assert(unique(pStartI(rStartI))==true)

%% 3. Construct design matrix (get event time histograms for each movement and task regressors 
% bumps for movement variables
[MOVE_base, ~, MOVE_func] = basis.log_cos(N_MOVE, RANGE_MOVE, DT, 10, false);
[TASK_base, ~, TASK_func] = basis.log_cos(N_TASK, RANGE_TASK, DT, 10, false);

% Design matrix (option 1): each trial type separately
for pt = 1:length(posTqTypes) % each position-torque pair
    tmpPosTqI = cell2mat(cellfun(@(a) strcmpi(a,posTqTypes{pt}), posTqC, 'un', 0 )) &pStartI'; % position - torque combination
    %tmpPosTqInull = ~cell2mat(cellfun(@(a) strcmpi(a,posTqTypes{pt}), posTqC, 'un', 0 )) &pStartI'; % tried concanenating both indices for a specific trials type and its rest, but such design could be problematic (overfitting), thus averted  
    dm_mv{pt,1} = 'move'; % regressor type 
    dm_mv{pt,2} = posTqTypes{pt}; % regressor name (position-torque combination) 
    %dm_mv{pt,3}(:,1) = histcounts([jkvt(tmpPosTqInull).rStartToPull],time_bin)'; % rest of the trials 
    dm_mv{pt,3} = histcounts([jkvt(tmpPosTqI).rStartToPull],time_bin)'; % discrete event histogram (consider '-1000 ms' to include movement preparatory period) 
    dm_mv{pt,4} = basis.conv(dm_mv{pt,3},MOVE_base); % convolution
end

% Design matrix (option 2): position and torque  
dm_mv1{1,1} = 'move'; 
dm_mv1{1,2} = 'position';
dm_mv1{1,3}(:,1) = histcounts([jkvt(leI).rStartToPull],time_bin)'; % left success trials
dm_mv1{1,3}(:,2) = histcounts([jkvt(riI).rStartToPull],time_bin)'; % right success trials
dm_mv1{1,4} = basis.conv(dm_mv1{1,3},MOVE_base); % convolution 

dm_mv1{2,1} = 'move'; 
dm_mv1{2,2} = 'torque';
dm_mv1{2,3}(:,1) = histcounts([jkvt(loI).rStartToPull],time_bin)'; % left success trials
dm_mv1{2,3}(:,2) = histcounts([jkvt(hiI).rStartToPull],time_bin)'; % right success trials
dm_mv1{2,4} = basis.conv(dm_mv1{2,3},MOVE_base); % convolution     

% task regressor: reward delivery 
dm_task{1,1} = 'task'; % regressor type 
dm_task{1,2} = 'rwd'; % regressor name
dm_task{1,3}(:,1) = histcounts([jkvt(~rwI).trEnd]+1000,time_bin)'; % not-rewarded
dm_task{1,3}(:,2) = histcounts([jkvt(rwI).rewardT],time_bin)'; % rewarded
dm_task{1,4} = basis.conv(dm_task{1,3},TASK_base); 

%% get bumps for continuous hand velocity variable and perform convolution
% generate bumps for reach and pull phases separately 
[~,k] = func.evtKinematics( jkvt, sessionDur );
hvcmS = k.handVel; % hand velocity in cm/S in ms
hvcmSb = hvcmS(1,time_bin(2:end)); %intm(hvcmS,length(time_bin)-1); % DO NOT USE 'intm' here! 

% reach phase
hvcmSbReach = max(0,hvcmSb); % reach phase (away from the initial position) 
hvcmSbReachBound = max(hvcmSbReach)*.9; % just use zero-to-max range
[hV_ReachBase, hV_ReachBins, hV_ReachFunc] = basis.linear_cos(N_HANDVEL, [0 hvcmSbReachBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to  
X_reachVel = hV_ReachFunc(hvcmSbReach); % convolution
% pull phase
absHvcmSbPull = abs(min(0,hvcmSb)); % pull phase (towards the initial position) 
absHvcmSbPullBound = max(absHvcmSbPull)*.9; % just use zero-to-absMax range
[hV_PullBase, hV_PullBins, hV_PullFunc] = basis.linear_cos(N_HANDVEL, [0 absHvcmSbPullBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to  
X_pullVel = hV_PullFunc(absHvcmSbPull); % convolution

dm_handVel{1,1} = 'continuous'; % regressor type
dm_handVel{1,2} = 'reachVel';   % regressor name
dm_handVel{1,3} = hvcmSbReach'; % continuous variable
dm_handVel{1,4} = X_reachVel; 

dm_handVel{2,1} = 'continuous'; % regressor type
dm_handVel{2,2} = 'pullVel';   % regressor name
dm_handVel{2,3} = absHvcmSbPull'; % continuous variable
dm_handVel{2,4} = X_pullVel; 

%% 4.2. Firing rate by task variables
%Direction-tuning: 79, 80, 35, 77, 85
i_cell = 103; % for loop

spike_time = clip(spkTimesCell{1,i_cell});
n_spike = length(spike_time);
[spike_bin,~,spike_binT] = histcounts(spike_time, time_bin);
spike_rate = sum(spike_bin) / (sessionDur/(1000)); % spikes/Sec

bin_size = 10; % ms
filter_sigma = 100; % ms
window = 5000; %5000; % ms
cut = 4 * filter_sigma / bin_size;

% align to the task event (reach start to pull)
[spike_rStart_count, time_tr_org] = alignToTaskEvent([jkvt(:).rStartToPull], spike_time, bin_size, bin_size, window + 4 * filter_sigma);      
spike_rStart = spike_rStart_count'.*(1000./binSize); % spike rate (Hz)
in_t = 1 + cut:length(time_tr_org) - cut;
time_tr = time_tr_org(in_t);

% binned spike counts aligned to trial end
[spike_rwdOrNo_count] = alignToTaskEvent([jkvt(:).trEnd]+1000, spike_time, bin_size, bin_size, window + 4 * filter_sigma);      
spike_rwdOrNo = spike_rwdOrNo_count'.*(1000./binSize); % spike rate (Hz)



spike_leLo = func.group_stat2(spike_rStart, leLoI(pStartI)+1); % mean PSTH, lelo vs. rest
spike_leHi = func.group_stat2(spike_rStart, leHiI(pStartI)+1); % mean PSTH, low vs. high 
spike_riLo = func.group_stat2(spike_rStart, riLoI(pStartI)+1); % mean PSTH, lelo vs. rest
spike_riHi = func.group_stat2(spike_rStart, riHiI(pStartI)+1); % mean PSTH, low vs. high 
spike_reward = func.group_stat2(spike_rwdOrNo, rwI+1); % mean PSTH, no Reward vs. reward

% normal filter
spike_leLo_conv = basis.normal_filter(spike_leLo, filter_sigma, bin_size);
spike_leHi_conv = basis.normal_filter(spike_leHi, filter_sigma, bin_size);
spike_riLo_conv = basis.normal_filter(spike_riLo, filter_sigma, bin_size);
spike_riHi_conv = basis.normal_filter(spike_riHi, filter_sigma, bin_size);
spike_reward_conv = basis.normal_filter(spike_reward, filter_sigma, bin_size); 

% log
leLo_log = log(spike_leLo_conv(in_t, :)/spike_rate);
leHi_log = log(spike_leHi_conv(in_t, :)/spike_rate);
riLo_log = log(spike_riLo_conv(in_t, :)/spike_rate);
riHi_log = log(spike_riHi_conv(in_t, :)/spike_rate);
reward_log = log(spike_reward_conv(in_t,:)/spike_rate); 

if PLOT
    figure(3); clf;
    subplot(3, 2, 1);
    plot(time_tr, leLo_log);
    title('leLo vs rest');
    xlabel('time from reachStart');
    axis tight
    
    subplot(3, 2, 2);
    plot(time_tr, leHi_log);
    title('leHi vs rest');
    xlabel('time from reachStart');
    axis tight

    subplot(3, 2, 3);
    plot(time_tr, riLo_log);
    title('riLo vs rest');
    xlabel('time from reachStart');
    axis tight

    subplot(3, 2, 4);
    plot(time_tr, riHi_log);
    title('riHi vs rest');
    xlabel('time from reachStart');
    axis tight
    
    subplot(3, 2, 5);
    plot(time_tr, reward_log);
    title('no reward vs reward');
    xlabel('time from trEnd');
    axis tight
end

%% 4. getting average response for parameter fitting (not necessary)
% 4.1. Firing rate by behavior
% coarse binning
%time_bin_s = (0:ceil(sessionDur/1000))';
%n_bin_s = length(time_bin_s);
ratio_time = floor(1 / (DT/1000)); % ratio of 1000ms to DT
ratio_time_100ms = floor(1 / (DT/100)); % ratio of 100ms to DT
coarse_100msBin = @(x) mean(reshape(x(1:floor(n_bin/ratio_time_100ms)*ratio_time_100ms), ratio_time_100ms, []))';   % mean with 100-ms bin
coarseMax_100msBin = @(x) max(reshape(x(1:floor(n_bin/ratio_time_100ms)*ratio_time_100ms), ratio_time_100ms, []))'; % max with 100-ms bin
spike_100ms = coarse_100msBin(spike_bin) * ratio_time; % spike rate 100-ms bin in Hz (convertion to Hz - multiply the ratio of 1000ms to the original binSize)
% bin max reach velocity
reachVel_100ms = coarseMax_100msBin(hvcmSbReach); 
reachVel_100msEdge = round(linspace(0,ceil(range(reachVel_100ms)*.9),11)); % to get 10 bins 
[rVel_mSpkR_100ms, rVel_sSpkR_100ms, rVel_bin_100ms] = func.group_stat(reachVel_100ms, spike_100ms, reachVel_100msEdge); % get the mean and sem spike rates per binned reach velocity 
rVel_logSpkRate_100ms = log(rVel_mSpkR_100ms / spike_rate); % 
% bin max pull velocity 
pullVel_100ms = coarseMax_100msBin(absHvcmSbPull); 
pullVel_100msEdge = round(linspace(0,ceil(range(pullVel_100ms)*.9),11)); 
[pVel_mSpkR_100ms, pVel_sSpkR_100ms, pVel_bin_100ms] = func.group_stat(pullVel_100ms, spike_100ms, pullVel_100msEdge); % get the mean and sem spike rates per binned reach velocity
pVel_logSpkRate_100ms = log(pVel_mSpkR_100ms / spike_rate); % 

if PLOT
    % plotting
    figure(1); clf;
    subplot(2, 1, 1);
    hold on;
    scatter(reachVel_100ms, spike_100ms, '.');
    errorbar(rVel_bin_100ms, rVel_mSpkR_100ms, rVel_sSpkR_100ms);
    subplot(2, 1, 2);
    plot(rVel_bin_100ms, rVel_logSpkRate_100ms);
end

if PLOT
    % plotting
    figure(2); clf;
    subplot(2, 1, 1);
    hold on;
    scatter(pullVel_100ms, spike_100ms, '.');
    errorbar(pVel_bin_100ms, pVel_mSpkR_100ms, pVel_sSpkR_100ms);
    subplot(2, 1, 2);
    plot(pVel_bin_100ms, pVel_logSpkRate_100ms);
end


%% 5. Parameter fitting
rStart_base = MOVE_func(time_tr); % base for movement covariates
reward_base = TASK_func(time_tr); % base for reward covariates
rVel_base = hV_ReachFunc(rVel_bin_100ms); 
pVel_base = hV_PullFunc(pVel_bin_100ms); 

% fitting weights by average response
rStart_proj = pinv(rStart_base' * rStart_base) * rStart_base';
reward_proj = pinv(reward_base' * reward_base) * reward_base';

w_leLo0 = rStart_proj * leLo_log; % left position, low torque
w_leHi0 = rStart_proj * leHi_log; % left position, high torque
w_riLo0 = rStart_proj * riLo_log; % right position, low torque
w_riHi0 = rStart_proj * riHi_log; % right position, high torque
w_reward0 = reward_proj * reward_log;
w_rVel0 = pinv(rVel_base' * rVel_base) * (rVel_base' * rVel_logSpkRate_100ms);
w_pVel0 = pinv(pVel_base' * pVel_base) * (pVel_base' * pVel_logSpkRate_100ms);
%w_h0 = pinv(h_base' * h_base) * (h_base' * spc_log); % spike history (auto correlation) is not included as the binsize = 10ms
w_c0 = log(spike_rate);

% rebuilding fitted kernel
leLo0 = rStart_base * w_leLo0;
leHi0 = rStart_base * w_leHi0;
riLo0 = rStart_base * w_riLo0;
riHi0 = rStart_base * w_riHi0;
reward0 = reward_base * w_reward0;
rVel0 = rVel_base * w_rVel0; % reach velocity
pVel0 = pVel_base * w_pVel0; % pull velocity

if PLOT
    figure(5); clf;
    subplot(4, 2, 1);
    hold on;
    plot(leLo_log); 
    plot(leLo0); 
    title('leLo vs rest');
    
    subplot(4, 2, 2);
    hold on;
    plot(leHi_log); 
    plot(leHi0); 
    title('leHi vs rest');
    
    subplot(4, 2, 3);
    hold on;
    plot(riLo_log); 
    plot(riLo0); 
    title('riLo vs rest');
    
    subplot(4, 2, 4);
    hold on;
    plot(riHi_log); 
    plot(riHi0); 
    title('riHi vs rest');
    
    subplot(4, 2, 5);
    hold on;
    plot(reward_log);
    plot(reward0);
    title('reward');   
    
    subplot(4, 2, 6); 
    hold on; 
    plot(rVel_bin_100ms, rVel_logSpkRate_100ms);
    plot(rVel_bin_100ms, rVel0);
    title('reach velocity')
    
    subplot(4, 2, 7); 
    hold on; 
    plot(pVel_bin_100ms, pVel_logSpkRate_100ms);
    plot(pVel_bin_100ms, pVel0);
    title('pull velocity')
end

%% 6. loss function
[X, prm] = func.normalize_add_constant([dm_mv(1,4)']); %, dm_task(1,4)']); %, dm_handVel(:,4)']); 

if CALC_PRM
    prm0 = [w_c0; w_leLo0(:)];% w_leHi0(:); w_riLo0(:); w_riHi0(:); w_reward0(:)]; % w_rVel0; w_pVel0];
else
    w_c0 = log(spike_rate);
    prm0 = [w_c0; rand(prm.n_var, 1) - 0.5];
end

prm0_norm = [prm0(1) + prm.mean * prm0(2:end); prm0(2:end) .* prm.std'];

lfunc = @(w) loss.log_poisson_loss(w, X, spike_bin', DT/1000); % make sure that DT here should be 1/1000, well not necessarily it depends on the binSize

%% optimization
algopts = {'algorithm','trust-region','Gradobj','on','Hessian','on', 'display', 'iter', 'maxiter', 100};
opts = optimset(algopts{:});
[prm1_norm, loss1, exitflag, output, grad, hessian] = fminunc(lfunc, prm0_norm, opts);

%% revert prm
prm1 = [prm1_norm(1) - (prm.mean ./ prm.std) * prm1_norm(2:end); prm1_norm(2:end) ./ prm.std']; % original
prm1_std_norm = sqrt(diag(inv(hessian)));
prm1_std = [prm1_std_norm(1); prm1_std_norm(2:end) ./ prm.std'];

%% plot leLo0 = rStart_base * w_leLo0;
leLo1 = rStart_base * reshape(prm1(prm.index{2}), [], 2);
leLo1_std = rStart_base * reshape(prm1_std(prm.index{2}), [], 2);

leHi1 = rStart_base * reshape(prm1(prm.index{3}), [], 2);
leHi1_std = rStart_base * reshape(prm1_std(prm.index{3}), [], 2);
 
riLo1 = rStart_base * reshape(prm1(prm.index{4}), [], 2);
riLo1_std = rStart_base * reshape(prm1_std(prm.index{2}), [], 2);

riHi1 = rStart_base * reshape(prm1(prm.index{5}), [], 2);
riHi1_std = rStart_base * reshape(prm1_std(prm.index{3}), [], 2);
 
reward1 = reward_base * reshape(prm1(prm.index{6}), [], 2);
reward1_std = reward_base * reshape(prm1_std(prm.index{4}), [], 2);
 
% rVel1 = rVel_base * prm1(prm.index{7});
% rVel1_std = rVel_base * prm1_std(prm.index{7});
% 
% pVel1 = pVel_base * prm1(prm.index{8});
% pVel1_std = pVel_base * prm1_std(prm.index{8});

%h1 = h_base * prm1(prm.index{3});
%h1_std = h_base * prm1_std(prm.index{3});

figure(6); clf;
subplot(4, 2, 1);
hold on;
plot(time_tr, leLo0);
plot(time_tr, leLo1);
title('leLo');

subplot(4, 2, 2);
hold on;
plot(time_tr, leHi0);
plot(time_tr, leHi1);
title('leHi');

subplot(4, 2, 3);
hold on;
plot(time_tr, riLo0);
plot(time_tr, riLo1);
title('riLo');

subplot(4, 2, 4);
hold on;
plot(time_tr, riHi0);
plot(time_tr, riHi1);
title('riHi');


%%
Xprm1 = X*prm1_norm; 
plot(Xprm1);
expXprm1 = exp(Xprm1); 
plot(expXprm1); 

%% individual unit psth aligned to a task event
% sort trials
tq = round([jkvt.pull_torque]./10); % torque 
tqT = tq(S.trI)'; tqT(:,2) = 1:length(S.trI); tqT(:,3) = S.trI;  sortByTq = sortrows(tqT,1); 
ps = [jkvt.reachP1]; % position 
psT = ps(S.trI)'; psT(:,2) = 1:length(S.trI); psT(:,3) = S.trI;  sortByPs = sortrows(psT,1); 
tqpsT(:,1) = tqT(:,1)+psT(:,1); tqpsT(:,2) = 1:length(S.trI); tqpsT(:,3) = S.trI; sortByTqPs = sortrows(tqpsT,1); % torque position combo

%load('binSpkCountSTRCTXWR40_082019.mat', 'rStartToPull')
%S = rStartToPull
i_cell =102; % 102 79, 80, 35, 77, 85
thisUnitSpkTimes = S.SpkTimes{i_cell}; 
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTq, i_cell, [3e3 2e3], [2e3 1.9e3]); 
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByPs, i_cell, [3e3 2e3], [2e3 1.9e3]);
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTqPs, i_cell, [3e3 2e3], [2e3 1.9e3]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [binSpkMatOut,binEdgesOut] = alignToTaskEvent(evt, spikeTime, binSize, stepSize, oneWindow)           
        binEdgesOut = -oneWindow:binSize:oneWindow; 
        bin1msEdge = cellfun(@(a) max(a-oneWindow,0):a+oneWindow+1, num2cell(evt),'un', 0); 
        bin1msCount = cellfun(@(a) histcounts(spikeTime,a), bin1msEdge, 'un',0); 
        binSpkMatOut = bin1msSpkCountMat(cell2mat(bin1msCount'),binSize,stepSize); 
end

function [unitIdSTC] = getUnitIdSpkTimesCell(spkTimesCellIn, sIn, unitIdS)
%unit Ids don't match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
% this function identifies the unit from 'S' in 'spkTimesCell' and outputs
% it as 'unitidSTC'

geomI = cell2mat(cellfun(@(a) isequal(a,sIn.geometry{unitIdS}), spkTimesCellIn(4,:), 'un', 0)); 
wfI = cell2mat(cellfun(@(a) isequal(a,sIn.meanWF{unitIdS}), spkTimesCellIn(6,:), 'un', 0)); 

unitIdSTC = find(geomI & wfI); 

end







