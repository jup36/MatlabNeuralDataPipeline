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
sessionDur = rwdTimes(end)+10*1000; % cut 10 seconds after the last reward
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
rwI = [jkvt(:).rewarded]'; 

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
dm_task{1,2} = 'rwd';  % regressor name
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
hvcmSbReachBound = max(hvcmSbReach); %*.8; % just use zero-to-max range
[hV_ReachBase, hV_ReachBins, hV_ReachFunc] = basis.linear_cos(N_HANDVEL, [0 hvcmSbReachBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to  
X_reachVel = hV_ReachFunc(hvcmSbReach); % convolution
% pull phase
absHvcmSbPull = abs(min(0,hvcmSb)); % pull phase (towards the initial position) 
absHvcmSbPullBound = max(absHvcmSbPull); %*.8; % just use zero-to-absMax range
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
[i_cell] = getUnitIdSpkTimesCell(spkTimesCell, S, 102); % get the unit index as in the 'spkTimesCell'
%i_cell = 102; % for loop

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
spike_rStart = spike_rStart_count'.*(1000./bin_size); % spike rate (Hz)
in_t = 1 + cut:length(time_tr_org) - cut;
time_tr = time_tr_org(in_t);

% binned spike counts aligned to trial end
[spike_rwdOrNo_count] = alignToTaskEvent([jkvt(:).trEnd]+1000, spike_time, bin_size, bin_size, window + 4 * filter_sigma);      
spike_rwdOrNo = spike_rwdOrNo_count'.*(1000./bin_size); % spike rate (Hz)

spike_leRi = func.group_stat2(spike_rStart, riI(pStartI)+1); % mean PSTH, lelo vs. rest
spike_loHi = func.group_stat2(spike_rStart, hiI(pStartI)+1); % mean PSTH, lelo vs. rest
spike_reward = func.group_stat2(spike_rwdOrNo, rwI+1); % mean PSTH, no Reward vs. reward

% normal filter
spike_leRi_conv = basis.normal_filter(spike_leRi, filter_sigma, bin_size);
spike_loHi_conv = basis.normal_filter(spike_loHi, filter_sigma, bin_size);
spike_reward_conv = basis.normal_filter(spike_reward, filter_sigma, bin_size); 

% log
leRi_log = log(spike_leRi_conv(in_t, :)/spike_rate);
loHi_log = log(spike_loHi_conv(in_t, :)/spike_rate);
reward_log = log(spike_reward_conv(in_t,:)/spike_rate); 

if PLOT
    figure(3); clf;
    subplot(2, 2, 1);
    plot(time_tr, leRi_log);
    title('left vs right');
    xlabel('time from reachStart');
    axis tight
    
    subplot(2, 2, 2);
    plot(time_tr, loHi_log);
    title('low vs high torque');
    xlabel('time from reachStart');
    axis tight
    
    subplot(2, 2, 3);
    plot(time_tr, reward_log);
    title('no reward vs reward');
    xlabel('time from trEnd');
    axis tight
end

%% 4. getting average response for parameter fitting (not necessary)
% 4.1. Firing rate by behavior
% coarse binning
ratio_time = floor(1 / (DT/1000)); % ratio of 1000ms to DT
ratio_time_100ms = floor(1 / (DT/100)); % ratio of 100ms to DT
coarse_100msBin = @(x) mean(reshape(x(1:floor(n_bin/ratio_time_100ms)*ratio_time_100ms), ratio_time_100ms, []))';   % mean with 100-ms bin
coarseMax_100msBin = @(x) max(reshape(x(1:floor(n_bin/ratio_time_100ms)*ratio_time_100ms), ratio_time_100ms, []))'; % max with 100-ms bin
spike_100ms = coarse_100msBin(spike_bin) * ratio_time; % spike rate 100-ms bin in Hz (convertion to Hz - multiply the ratio of 1000ms to the original binSize)
%coarse_Bin = @(x) mean(reshape(x(1:floor(n_bin/ratio_time)*ratio_time), ratio_time, []))';   % mean with 1000-ms bin
%coarseMax_Bin = @(x) max(reshape(x(1:floor(n_bin/ratio_time)*ratio_time), ratio_time, []))'; % max with 1000-ms bin
%spike_1000ms = coarse_Bin(spike_bin) * ratio_time; % spike rate 1000-ms bin in Hz (convertion to Hz - multiply the ratio of 1000ms to the original binSize)
% bin max reach velocity
reachVel_100ms = coarseMax_100msBin(hvcmSbReach); 
reachVel_100msEdge = round(linspace(0,ceil(range(reachVel_100ms)),11)); % to get 10 bins %round(linspace(0,ceil(range(reachVel_100ms)*.8),11)); % to get 10 bins 
[rVel_mSpkR_100ms, rVel_sSpkR_100ms, rVel_bin_100ms] = func.group_stat(reachVel_100ms, spike_100ms, reachVel_100msEdge); % get the mean and sem spike rates per binned reach velocity 
rVel_logSpkRate_100ms = log(rVel_mSpkR_100ms / spike_rate); % 
%reachVel = coarseMax_Bin(hvcmSbReach); 
%reachVel_Edge = round(linspace(0,ceil(range(reachVel)*.9),6)); % to get 5 bins 
%[rVel_mSpkR, rVel_sSpkR, rVel_bin] = func.group_stat(reachVel, spike_1000ms, reachVel_Edge); % get the mean and sem spike rates per binned reach velocity 
%rVel_logSpkRate = log(rVel_mSpkR / spike_rate); % 

% bin max pull velocity 
pullVel_100ms = coarseMax_100msBin(absHvcmSbPull); 
pullVel_100msEdge = round(linspace(0,ceil(range(pullVel_100ms)),11)); %round(linspace(0,ceil(range(pullVel_100ms)*.8),11)); 
[pVel_mSpkR_100ms, pVel_sSpkR_100ms, pVel_bin_100ms] = func.group_stat(pullVel_100ms, spike_100ms, pullVel_100msEdge); % get the mean and sem spike rates per binned reach velocity
pVel_logSpkRate_100ms = log(pVel_mSpkR_100ms / spike_rate); % 
%pullVel = coarseMax_Bin(absHvcmSbPull); 
%pullVel_Edge = round(linspace(0,ceil(range(pullVel)*.9),11)); 
%[pVel_mSpkR, pVel_sSpkR, pVel_bin] = func.group_stat(pullVel, spike_1000ms, pullVel_Edge); % get the mean and sem spike rates per binned reach velocity
%pVel_logSpkRate = log(pVel_mSpkR / spike_rate); % 

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

w_leRi0 = rStart_proj * leRi_log; % w0 for the position (left vs right) variable
w_loHi0 = rStart_proj * loHi_log; % w0 for the torque (low vs high) variable
w_reward0 = reward_proj * reward_log;
w_rVel0 = pinv(rVel_base' * rVel_base) * (rVel_base' * rVel_logSpkRate_100ms);
w_pVel0 = pinv(pVel_base' * pVel_base) * (pVel_base' * pVel_logSpkRate_100ms);
%w_h0 = pinv(h_base' * h_base) * (h_base' * spc_log); % spike history (auto correlation) is not included as the binsize = 10ms
w_c0 = log(spike_rate);

% rebuilding fitted kernel
leRi0 = rStart_base * w_leRi0;
loHi0 = rStart_base * w_loHi0;
reward0 = reward_base * w_reward0;
rVel0 = rVel_base * w_rVel0; % reach velocity
pVel0 = pVel_base * w_pVel0; % pull velocity

if PLOT
    figure(5); clf;
    subplot(3, 2, 1);
    hold on;
    plot(leRi_log); 
    plot(leRi0); 
    title('left vs right');
    
    subplot(3, 2, 2);
    hold on;
    plot(loHi_log); 
    plot(loHi0); 
    title('low vs high');
    
    subplot(3, 2, 3);
    hold on;
    plot(reward_log);
    plot(reward0);
    title('no-reward vs reward');   
    
    subplot(3, 2, 4); 
    hold on; 
    plot(rVel_bin_100ms, rVel_logSpkRate_100ms);
    plot(rVel_bin_100ms, rVel0);
    title('reach velocity')
    
    subplot(3, 2, 5); 
    hold on; 
    plot(pVel_bin_100ms, pVel_logSpkRate_100ms);
    plot(pVel_bin_100ms, pVel0);
    title('pull velocity')
end

%% 6. loss function
[X, prm] = func.normalize_add_constant([dm_mv1(:,4)', dm_task(1,4)', dm_handVel(:,4)']); 

if CALC_PRM
    prm0 = [w_c0; w_leRi0(:); w_loHi0(:); w_reward0(:); w_rVel0; w_pVel0];
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
leRi1 = rStart_base * reshape(prm1(prm.index{2}), [], 2);
leRi1_std = rStart_base * reshape(prm1_std(prm.index{2}), [], 2);

loHi1 = rStart_base * reshape(prm1(prm.index{3}), [], 2);
loHi1_std = rStart_base * reshape(prm1_std(prm.index{3}), [], 2);
 
reward1 = reward_base * reshape(prm1(prm.index{4}), [], 2);
reward1_std = reward_base * reshape(prm1_std(prm.index{4}), [], 2);
 
rVel1 = rVel_base * prm1(prm.index{5});
rVel1_std = rVel_base * prm1_std(prm.index{5});
 
pVel1 = pVel_base * prm1(prm.index{6});
pVel1_std = pVel_base * prm1_std(prm.index{6});

figure(6); clf;
subplot(3, 2, 1);
hold on;
plot(time_tr, leRi0);
plot(time_tr, leRi1);
title('left vs right');

subplot(3, 2, 2);
hold on;
plot(time_tr, loHi0);
plot(time_tr, loHi1);
title('low vs high');

subplot(3, 2, 3);
hold on;
plot(time_tr, reward0);
plot(time_tr, reward1);
title('no-reward vs reward');

subplot(3, 2, 4);
hold on;
plot(rVel_bin_100ms, rVel0);
plot(rVel_bin_100ms, rVel1);
title('reach velocity');

subplot(3, 2, 5);
hold on;
plot(pVel_bin_100ms, pVel0);
plot(pVel_bin_100ms, pVel1);
title('pull velocity');

%% get model prediction firing rates
Xprm1 = X*prm1_norm; 
%figure; plot(Xprm1);
expXprm1 = basis.normal_filter(exp(Xprm1), filter_sigma, 10); 

spike_bin_conv = basis.normal_filter((spike_bin.*(1000/DT))', filter_sigma, 10); % filter_sigma = 100 ms
figure; hold on; plot(spike_bin_conv); plot(expXprm1); hold off;  

% get model prediction PSTHs
model_rStart = func.alignGlmOutToTaskEvent([jkvt(:).rStartToPull], time_bin, expXprm1, 10, window + 4 * filter_sigma); % fitted rate aligned to reachStart  
model_rwdOrNo = func.alignGlmOutToTaskEvent([jkvt(:).trEnd]+1000, time_bin, expXprm1, 10, window + 4 * filter_sigma);  % fitted rate aligned to trEnd+1000ms

model_leRi = func.group_stat2(model_rStart, riI(pStartI)+1); % fitted mean left vs. mean right
model_loHi = func.group_stat2(model_rStart, hiI(pStartI)+1); % fitted mean low vs. mean high
model_reward = func.group_stat2(model_rwdOrNo, rwI+1); % fitted mean no-reward vs. mean reward

%% Variance analysis (R^2)
r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
r2_entire = r2(spike_bin_conv, expXprm1); 
% calculate only within PSTHs? Not across the whole session
model_rStart_psth2s = func.alignGlmOutToTaskEvent([jkvt(:).rStartToPull], time_bin, expXprm1, 10, 2000); % fitted rate aligned to reachStart  
spike_rStart_psth2s = func.alignGlmOutToTaskEvent([jkvt(:).rStartToPull], time_bin, spike_bin_conv, 10, 2000); % fitted rate aligned to reachStart  
r2_psth = r2(reshape(spike_rStart_psth2s,[],1),reshape(model_rStart_psth2s,[],1));

%% decoding target position
% p(r|x=left,x~)
% align fitted data to task events
% posterior probability (target position, left)
% take the spike train of the trial
aSpkTrain = basis.normal_filter(spike_rStart(:,78), filter_sigma, bin_size);
figure; plot(aSpkTrain)

postLe = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(aSpkTrain), model_leRi(:,1), 'un', 0)); 
postRi = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(aSpkTrain), model_leRi(:,2), 'un', 0)); 

%figure; plot(postLe./(postLe+postRi)); 
%figure; plot(model_leRi_conv); 

% posterior probability (torque, low)
postLo = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(aSpkTrain), model_loHi(:,1), 'un', 0)); 
postHi = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(aSpkTrain), model_loHi(:,2), 'un', 0)); 

figure; plot(postLo./(postLo+postHi)); 
%figure; plot(model_loHi_conv); 



%% individual unit psth aligned to a task event
% sort trials
tq = round([jkvt.pull_torque]./10); % torque 
tqT = tq(S.trI)'; tqT(:,2) = 1:length(S.trI); tqT(:,3) = S.trI;  sortByTq = sortrows(tqT,1); 
ps = [jkvt.reachP1]; % position 
psT = ps(S.trI)'; psT(:,2) = 1:length(S.trI); psT(:,3) = S.trI;  sortByPs = sortrows(psT,1); 
tqpsT(:,1) = tqT(:,1)+psT(:,1); tqpsT(:,2) = 1:length(S.trI); tqpsT(:,3) = S.trI; sortByTqPs = sortrows(tqpsT,1); % torque position combo

%load('binSpkCountSTRCTXWR40_082019.mat', 'rStartToPull')
%S = rStartToPull
cellI =102; % 102 79, 80, 35, 77, 85
thisUnitSpkTimes = S.SpkTimes{cellI}; 
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTq, cellI, [3e3 2e3], [2e3 1.9e3]); 
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByPs, cellI, [3e3 2e3], [2e3 1.9e3]);
individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTqPs, cellI, [3e3 2e3], [2e3 1.9e3]);

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

function [winModelRate] = alignGlmOutToTaskEvent(evt, time_bin, modelRate, binSize, oneWindow)
    windowBins = [floor(-oneWindow/binSize)+1,ceil(oneWindow/binSize)]; 
    binEvt = find(histcounts(evt,time_bin)'); % identify bins per event
    windowBound = arrayfun(@(a) a+windowBins,binEvt,'un',0); 
    for t = 1:length(windowBound)
        if windowBound{t}(1)>=1 && windowBound{t}(2)<=length(modelRate)
            winModelRateC{t} = modelRate(windowBound{t}(1):windowBound{t}(2)); 
        else
            winModelRateC{t} = NaN(sum(abs(windowBins))+1,1); 
        end
    end
    %winModelRateC = cellfun(@(a) modelRate(a(1):a(2)), windowBound, 'un',0); 
    winModelRate = cell2mat(winModelRateC);   
end







