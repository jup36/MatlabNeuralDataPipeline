
clc; clearvars; close all;

%% 1. Load data
filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'; 
cd(filePath)
% neural data
% neural and behavioral data
spkDir = dir('binSpkCountSTRCTX*');
load(fullfile(spkDir(1).folder, spkDir(1).name),'spkTimesCell','jkvt','meta')
%S=rStartToPull; 

% task parameter
DT = 20; % 20 ms (width of timeBin)

N_MOVE = 10; % reach start, n bumps
RANGE_MOVE = [0 3000]; % ms (-1000 to 2000 ms relative to rStartToPull)
BIAS_MOVE = 5;

N_TASK = 10; % reward
RANGE_TASK = 3000; % ms (0 to 3000 ms relative to reward delivery) 
BIAS_TASK = 8;

N_H = 10; % n bumps for spike history
RANGE_H = 200; 
BIAS_H = 10;

% code parameter
PLOT = true;
CALC_PRM = false;

%% 2. Behavioral data 
% detect events and continuous joystick pull trajectory  
sessionDur = round(str2double(meta.fileTimeSecs)*1000); % in ms 
session_range = [1 sessionDur]; 
clip = @(t) t(session_range(1) <= t & t <= session_range(2)) - session_range(1);
time_bin = (0:DT:sessionDur)';
n_bin = length(time_bin) - 1;

rStartI = cellfun(@(c) ~isempty(c), {jkvt(:).rStartToPull});
pullStartI = cellfun(@(c) ~isempty(c), {jkvt(:).pullStarts}); 
rewardTimeI = cellfun(@(c) ~isnan(c), {jkvt(:).rewardT}); 

% if reach Start is missing, put estimated values based on the pullStarts with (-200 ms offset)
[jkvt(~rStartI&pullStartI).rStartToPull] = deal([jkvt(~rStartI&pullStartI).pullStarts]-200); 
rStartI(~rStartI&pullStartI) = true; 
assert(unique(pullStartI(rStartI))==true)

%% get event time histograms for each movement and task regressors 
% parse joystick position - torque combinations
[posTqC,posTqTypes] = posTrqCombinations( jkvt ); 
for pt = 1:length(posTqTypes)
    tmpPosTqI = cell2mat(cellfun(@(a) strcmpi(a,posTqTypes{pt}), posTqC, 'un', 0 )); % position - torque combination
    dm_move{pt,1} = 'move'; % regressor type 
    dm_move{pt,2} = posTqTypes{pt}; % regressor name (position-torque combination) 
    dm_move{pt,3} = histcounts([jkvt(tmpPosTqI).rStartToPull]-1000,time_bin)';  % discrete event histogram, '-1000 ms' to include movement preparatory period 
end

% task regressor: joystick ready (trial start) 
dm_task{1,1} = 'task'; % regressor type
dm_task{1,2} = 'jsReady'; % regressor name
dm_task{1,3} = histcounts([jkvt(:).trJsReady],time_bin)'; % discrete event histogram

% task regressor: reward delivery 
dm_task{2,1} = 'task'; % regressor type 
dm_task{2,2} = 'rwd'; % regressor name
dm_task{2,3} = histcounts([jkvt(:).rewardT],time_bin); % discrete event histogram

%% get bumps (timed consine functions) for discrete MOVE and TASK variables using linear or log cosine functions and perform convolution 
% bumps (linear consine) for movement variables (aligned to -1000ms relative to reachStarts to include preparatory period)
[MOVE_base, ~, ~] = basis.linear_cos(N_MOVE, RANGE_MOVE, DT, false);
% convolution of each movement variable 
for t = 1:size(dm_move,1)
    dm_move{t,4} = basis.conv(dm_move{t,3},MOVE_base); 
end
clearvars t

[TASK_base, ~, ~] = basis.log_cos(N_TASK, RANGE_TASK, DT, 10, false);
% convolution of each task variable
for t = 1:size(dm_task)
    dm_task{t,4} = basis.conv(dm_task{t,3},TASK_base); 
end
clearvars t

%% get bumps for continuous hand velocity variable and perform convolution
% generate bumps for reach and pull phases separately 
[~,k] = func.evtKinematics( jkvt, sessionDur );
hvcmS = k.handVel; % hand velocity in cm/S
hvcmSb = intm(hvcmS,length(time_bin)-1); % bin (skip smoothing)
% reach phase
hvcmSbReach = max(0,hvcmSb); % reach phase (away from the initial position) 
hvcmSbReachBound = max(hvcmSbReach); % just use zero-to-max range
[hV_ReachBase, hV_ReachBins, hV_ReachFunc] = basis.linear_cos(N_HANDVEL, [0 hvcmSbReachBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to  
X_reachVel = hV_ReachFunc(hvcmSbReach); % convolution
% pull phase
absHvcmSbPull = abs(min(0,hvcmSb)); % pull phase (towards the initial position) 
absHvcmSbPullBound = max(absHvcmSbPull); % just use zero-to-absMax range
[hV_PullBase, hV_PullBins, hV_PullFunc] = basis.linear_cos(N_HANDVEL, [0 absHvcmSbPullBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to  
X_pullVel = hV_ReachFunc(absHvcmSbPull); % convolution

dm_handVel{1,1} = 'continuous'; % regressor type
dm_handVel{1,2} = 'handVel'; % regressor name
dm_handVel{1,3} = hvcmSb; % continuous variable

%% 3. Spike
% spike
%Direction-tuning: 79, 80, 35, 77, 85
i_cell = 79; % for loop

spike_time = clip(spkTimesCell{1,i_cell});
n_spike = length(spike_time);
[spike_bin,~,spike_binT] = histcounts(spike_time, time_bin);
spike_bin_dT = diff([0; spike_binT]); 
spike_bin(spike_binT(spike_bin_dT<=2))=0; 
spike_bin = min(spike_bin,1); 
spike_rate = sum(spike_bin) / (sessionDur/1000); % spikes/Sec

% spike bump
[h_base, h_time] = basis.log_cos(N_H, [DT, RANGE_H], DT, 10, true);
%[h_base, h_time] = basis.log_cos(n_h_bump, [1, n_h_time], 1, 10);
%[h_base, h_time] = basis.log_cos(N_H, [DT, RANGE_H], DT, BIAS_H, false); %true);

% convolution
X_h = basis.conv(min(spike_bin,1), h_base, h_time(1)/DT);

%% 4. getting average response for parameter fitting (not necessary)
% 4.1. Firing rate by behavior
% coarse binning
time_bin_s = (0:ceil(sessionDur/1000))';
n_bin_s = length(time_bin_s);
ratio_time = floor(1 / (1/1000));
coarse_bin = @(x) mean(reshape(x(1:floor(n_bin/ratio_time)*ratio_time), ratio_time, []))'; % spike rate 1-s bin
coarseMax_bin = @(x) max(reshape(x(1:floor(n_bin/ratio_time)*ratio_time), ratio_time, []))'; % max binning with 1-s bin

spike_s = coarse_bin(spike_bin) * ratio_time;
speed_s = coarseMax_bin(jsVcmS); 
speed_edge = 0:ceil(range(speed_s)); 
[speed_mean, speed_sem, speed_bin] = func.group_stat(speed_s, spike_s, speed_edge);
speed_log = log(speed_mean / spike_rate);

work_s = coarseMax_bin(jsWuJ);
work_edge = 0:round(range(work_s)/11):ceil(range(work_s)); 
[work_mean, work_sem, work_bin] = func.group_stat(work_s, spike_s, work_edge);
work_log = log(work_mean / spike_rate);

if PLOT
    % plotting
    figure(1); clf;
    subplot(2, 1, 1);
    hold on;
    scatter(speed_s, spike_s, '.');
    errorbar(speed_bin, speed_mean, speed_sem);
    subplot(2, 1, 2);
    plot(speed_bin, speed_log);
end

if PLOT
    % plotting
    figure(2); clf;
    subplot(2, 1, 1);
    hold on;
    scatter(work_s, spike_s, '.');
    errorbar(work_bin, work_mean, work_sem);
    subplot(2, 1, 2);
    plot(work_bin, work_log);
end

%% 4.2. Firing rate by task variables
bin_size = 10; % ms
filter_sigma = 100; % ms
window = 5000; % ms
cut = 4 * filter_sigma / bin_size;

% binned spike counts aligned to trial start 
[spike_rStart_ms, time_tr] = func.fast_align([jkvt(:).rStartToPull], spike_time, ...
    bin_size, window + 4 * filter_sigma); %  ms to include preparatory period spike activity
in_t = 1 + cut:length(time_tr) - cut;
time_tr = time_tr(in_t);
spike_rStart = spike_rStart_ms.*ratio_time; % multiply by ratio_time 1000 to convert to Hz

% binned spike counts aligned to trial end
spike_trEnd_ms = func.fast_align([jkvt(:).trEnd], spike_time, ...
    bin_size, window + 4 * filter_sigma);
spike_trEnd = spike_trEnd_ms.*ratio_time; % multiply by ratio_time 1000 to convert to Hz

% mean PSTH trials by direction, by torque
leRiI = zeros(size(jkvt,2), 1); 
leRiI(pullStartI & leftI) = 1; % left
leRiI(pullStartI & rightI) = 2; % right

loHiI = zeros(size(jkvt,2), 1);
loHiI(pullStartI & lowTqI) = 1; % low torque
loHiI(pullStartI & highTqI) = 2; % high torque

spike_leRi = func.group_stat2(spike_rStart, leRiI(pullStartI)); % mean PSTH, left vs. right
spike_loHi = func.group_stat2(spike_rStart, loHiI(pullStartI)); % mean PSTH, low vs. high 
spike_reward = func.group_stat2(spike_trEnd, rwI+1); % mean PSTH, no Reward vs. reward

% normal filter
spike_leRi_conv = basis.normal_filter(spike_leRi, filter_sigma, bin_size);
spike_loHi_conv = basis.normal_filter(spike_loHi, filter_sigma, bin_size);
spike_reward_conv = basis.normal_filter(spike_reward, filter_sigma, bin_size);

leRi_log = log(spike_leRi_conv(in_t, :)/spike_rate);
loHi_log = log(spike_loHi_conv(in_t, :)/spike_rate); 
reward_log = log(spike_reward_conv(in_t,:)/spike_rate); 

if PLOT
    figure(3); clf;
    subplot(3, 1, 1);
    plot(time_tr, leRi_log);
    title('left vs right');
    xlabel('time from reachStart');
    axis tight

    subplot(3, 1, 2);
    plot(time_tr, loHi_log);
    title('low tq vs high tq');
    xlabel('time from reachStart');
    axis tight
    
    subplot(3, 1, 3);
    plot(time_tr, reward_log);
    title('no reward vs reward');
    xlabel('time from trEnd');
    axis tight
end

%% 4.3. Average autocorrelogram
% get autocorrelation
%[spc, time_spc] = xcorr(spike_bin,200);
%[spc, time_spc] = func.cross_corr(spike_time, spike_time, DT, RANGE_H); % cross_corr.c returns crosscorr values normalized in a way that the 0 lag value to be 1, iow, crosscorr relative to the total spike counts
[spc, time_spc_s] = func.cross_corr(spike_time/1000, spike_time/1000, 1/1000, RANGE_H(end)/1000);
time_spc_s = time_spc_s(RANGE_H(end)/DT+2:end);
time_spc_ms = round(time_spc_s*1000); 
spc = spc((RANGE_H(end)/DT)+2:end); % take one side
spc_log = log(spc / spike_rate + exp(-10));

if PLOT
    figure(4); clf;
    subplot(2, 1, 1);
    hold on;
    plot(time_spc_ms, spc);
    plot(time_spc_ms([1, end]), [spike_rate, spike_rate], 'k:');

    subplot(2, 1, 2);
    plot(time_spc_ms, spc_log);
end

%% 5. Parameter fitting
start_base_w = start_func(time_tr);
reward_base = reward_func(time_tr);
speed_base = speed_func(speed_bin);

% fitting weights by average response
task_proj = pinv(start_base_w' * start_base_w) * start_base_w';
w_leRi0 = task_proj * leRi_log;
w_loHi0 = task_proj * loHi_log;
w_reward0 = task_proj * reward_log;
w_speed0 = pinv(speed_base' * speed_base) * (speed_base' * speed_log);
w_h0 = pinv(h_base' * h_base) * (h_base' * spc_log);
w_c0 = log(spike_rate);

% rebuilding fitted kernel
leRi0 = start_base_w * w_leRi0;
loHi0 = start_base_w * w_loHi0;
reward0 = start_base_w * w_reward0;
speed0 = speed_base * w_speed0;
h0 = h_base * w_h0;

if PLOT
    figure(5); clf;
    subplot(3, 2, 1);
    hold on;
    plot(speed_bin, speed_log);
    plot(speed_bin, speed0);
    title('speed');
    
    subplot(3, 2, 2);
    hold on;
    plot(spc_log);
    plot(h0)
    title('spike history');
    
    subplot(3, 2, 3);
    hold on;
    plot(leRi_log);
    plot(leRi0);
    title('leRi Dir');
    
    subplot(3, 2, 4);
    hold on;
    plot(loHi_log);
    plot(loHi0);
    title('loHi tq');
    
    subplot(3, 2, 5);
    hold on;
    plot(reward_log);
    plot(reward0);
    title('reward');   
end

%% 6. loss function
%[X, prm] = func.normalize_add_constant({X_leRi, X_loHi, X_reward, X_speed, X_h}, false);
%prm0 = [w_c0; w_leRi0(:); w_loHi0(:); w_reward0(:); w_speed0; w_h0];
[X, prm] = func.normalize_add_constant({X_leRi, X_h}, false);
%[X, prm] = func.normalize_add_constant({X_h}, false);

if CALC_PRM
    prm0 = [w_c0; w_leRi0(:); w_loHi0(:); w_reward0(:); w_speed0];
else
    w_c0 = log(spike_rate);
    prm0 = [w_c0; rand(prm.n_var, 1) - 0.5];
end

prm0_norm = [prm0(1) + prm.mean * prm0(2:end); prm0(2:end) .* prm.std'];

lfunc = @(w) loss.log_poisson_loss(w, X, spike_bin, .001); % make sure that DT here should be 1/1000

%% optimization
algopts = {'algorithm','trust-region','Gradobj','on','Hessian','on', 'display', 'iter', 'maxiter', 100};
opts = optimset(algopts{:});
[prm1_norm, loss1, exitflag, output, grad, hessian] = fminunc(lfunc, prm0_norm, opts);

%% revert prm
prm1 = [prm1_norm(1) - (prm.mean ./ prm.std) * prm1_norm(2:end); prm1_norm(2:end) ./ prm.std'];
prm1_std_norm = sqrt(diag(inv(hessian)));
prm1_std = [prm1_std_norm(1); prm1_std_norm(2:end) ./ prm.std'];

%% plot
leRi1 = start_base_w * reshape(prm1(prm.index{2}), [], 2);
leRi1_std = start_base_w * reshape(prm1_std(prm.index{2}), [], 2);

loHi1 = start_base_w * reshape(prm1(prm.index{3}), [], 2);
loHi1_std = start_base_w * reshape(prm1_std(prm.index{3}), [], 2);
 
reward1 = reward_base * reshape(prm1(prm.index{4}), [], 2);
reward1_std = reward_base * reshape(prm1_std(prm.index{4}), [], 2);
 
speed1 = speed_base * prm1(prm.index{5});
speed1_std = speed_base * prm1_std(prm.index{5});

h1 = h_base * prm1(prm.index{3});
h1_std = h_base * prm1_std(prm.index{3});

figure(6); clf;
subplot(3, 2, 1);
hold on;
plot(time_tr, leRi0);
plot(time_tr, leRi1);
%for i = 1:2
%    errorbar(time_tr(1:10:end), leRi1(1:10:end, i), leRi1_std(1:10:end, i));
%end
title('Dir leRi');

subplot(3, 2, 2);
hold on;
plot(time_tr, loHi0);
plot(time_tr, loHi1);
% for i_task = 1:2
%     errorbar(time_task(1:10:end), choice1(1:10:end, i_task), choice1_std(1:10:end, i_task));
% end
title('Tq loHi');

subplot(3, 2, 3);
hold on;
plot(speed_bin, speed0);
plot(speed_bin, speed1);
title('speed');

subplot(3, 2, 4);
hold on;
plot(time_tr, reward0);
plot(time_tr, reward1);
% for i = 1:2
%     errorbar(time_tr(1:10:end), reward1(1:10:end,i), reward1_std(1:10:end,i));
% end
title('reward');

subplot(3, 2, 5);
hold on;
plot(h_time, h0);
errorbar(h_time(1:end), h1(1:end), h1_std(1:end));
title('spike history');

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

%unit = 80; 
thisUnitSpkTimes = S.SpkTimes{i_cell}; 
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTq, i_cell, [3e3 2e3], [2e3 1.9e3]); 
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByPs, i_cell, [3e3 2e3], [2e3 1.9e3]);
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTqPs, i_cell, [3e3 2e3], [2e3 1.9e3]);











