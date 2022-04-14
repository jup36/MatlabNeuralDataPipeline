%% 'glm_reachPull_leRi_loHi_cvR2.m' fits glm (linear poisson model) with 5-fold cross-validation, 
% which fits the neuronal spike rate with a set of discrete and continuous
% movement and task variables. 
% Here's the list of variables in the model: 
% Constant regressor: 
%  constant
% Discrete regressors: 
%  Target position: left vs. right reachStarts (range: -2 to 2s relative to reachStarts)
%  Joystick load:   low vs. high load reachStarts (range: -2 to 2s)
%  Reward:          no-reward vs. reward (range: 0 to 4s)
% Continuous regressors: 
%  Reach speed: a set of 10 regressors corresponding to 10 reach speed intervals
%  Pull speed:  a set of 10 regressors corresponding to 10 pull speed intervals
% 5-fold cross-validation ('crossValTrainTestSetsEachTrialType.m'): 
%  The model is fitted on the four folds of data, and is evaluated with the rest of the data. 
%  For optimization, fractions of data corresponding to the test set are excluded ('excludeSegments.m').  
% Evaluation: 
%  Goodness-of-fit is evaluated on the test set.
%  To estimate the upper bound of explanatory power of each regressor, R^2
%  is computed for the model fitted with each regressor alone. 
%  To estimate the lower bound, R^2 is computed for the model fitted with
%  excluding each regressor only. 

%% 0. Load data
cd(filePath)
spkDir = dir('binSpkCountSTRCTX*');
load(fullfile(spkDir(1).folder, spkDir(1).name),'spkTimesCell','jkvt')
%S=rStartToPull; 

%% 1. Task parameters
% task parameter
p.DT = 10; % 10 ms (width of timeBin)
p.N_MOVE = 20; % reach start, n bumps
p.RANGE_MOVE = [0 4000];
p.BIAS_MOVE = 5;
p.N_TASK = 20; % reward
p.RANGE_TASK = [0 4000]; % ms (0 to 3000 ms relative to reward delivery)
p.BIAS_TASK = 8;
p.N_HANDVEL = 10;
p.filter_sigma = 100; % ms
p.window = 5000; %5000; % ms
p.cut = 4 * p.filter_sigma / p.DT;
r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
p.isStr = cell2mat(spkTimesCell(5,:));

% code parameter
PLOT = false;
CALC_PRM = true;

%% 2. Behavioral data
% detect events and continuous joystick pull trajectory
%rwdTimes = cell2mat(cellfun(@(a) a(~isnan(a)), {jkvt(:).rewardT}, 'un', 0));
trEnds = cell2mat(cellfun(@(a) a(~isnan(a)), {jkvt(:).trEnd}, 'un', 0));
sessionDur = trEnds(end)+10*1000; % cut 10 seconds after the last reward
session_range = [1 sessionDur];
clip = @(t) t(session_range(1) <= t & t <= session_range(2)) - session_range(1);
time_bin = (0:p.DT:sessionDur)';
n_bin = length(time_bin) - 1;

% trial indices
[posTqC,posTqTypes] = posTrqCombinations( jkvt );
rStartI = cellfun(@(c) ~isempty(c), {jkvt(:).rStartToPull})';
pStartI = cellfun(@(c) ~isempty(c), {jkvt(:).pullStarts})';
leI = cell2mat(cellfun(@(a) contains(a,'p1'), posTqC, 'un', 0 )) & pStartI;
riI = cell2mat(cellfun(@(a) contains(a,'p2'), posTqC, 'un', 0 )) & pStartI;
loI = cell2mat(cellfun(@(a) contains(a,'t1'), posTqC, 'un', 0 )) & pStartI;
hiI = cell2mat(cellfun(@(a) contains(a,'t2'), posTqC, 'un', 0 )) & pStartI;
rwI = [jkvt(:).rewarded]';

% if reach Start is missing, put estimated values based on the pullStarts with (-200 ms offset)
if sum(~rStartI & pStartI)>=1
    tempMissingRstarts = find(~rStartI & pStartI);
    for tt = 1:length(tempMissingRstarts)
        jkvt(tempMissingRstarts(tt)).rStartToPull = jkvt(tempMissingRstarts(tt)).pullStarts - 200;
    end
end
rStartI = cellfun(@(c) ~isempty(c), {jkvt(:).rStartToPull})';
assert(unique(pStartI(rStartI))==true)

% To save task variables
tV.sessionDur = sessionDur;
tV.clip = clip;
tV.time_bin = time_bin;
tV.evtRstart = [jkvt(:).rStartToPull];
tV.evtRwd = [jkvt(:).trEnd]+1000;
tV.posTqC = posTqC;
tV.posTqTypes = posTqTypes;
tV.rStartI = rStartI;
tV.pStartI = pStartI;
tV.leI = leI;
tV.riI = riI;
tV.loI = loI;
tV.hiI = hiI;
tV.rwI = rwI;

% bumps for movement variables
[MOVE_base, ~, MOVE_func] = basis.log_cos(p.N_MOVE, p.RANGE_MOVE, p.DT, 10, false);
[TASK_base, ~, TASK_func] = basis.log_cos(p.N_TASK, p.RANGE_TASK, p.DT, 10, false);

% Discrete regressors position (full)
dmFull_mv1{1,1} = 'move';
dmFull_mv1{1,2} = 'position';
dmFull_mv1{1,3}(:,1) = histcounts([jkvt(leI).rStartToPull]-2000,time_bin)'; % left success trials
dmFull_mv1{1,3}(:,2) = histcounts([jkvt(riI).rStartToPull]-2000,time_bin)'; % right success trials
dmFull_mv1{1,4} = basis.conv(dmFull_mv1{1,3},MOVE_base); % convolution
% Discrete regressors torque (full)
dmFull_mv1{2,1} = 'move';
dmFull_mv1{2,2} = 'torque';
dmFull_mv1{2,3}(:,1) = histcounts([jkvt(loI).rStartToPull]-2000,time_bin)'; % low torque success trials
dmFull_mv1{2,3}(:,2) = histcounts([jkvt(hiI).rStartToPull]-2000,time_bin)'; % high torque success trials
dmFull_mv1{2,4} = basis.conv(dmFull_mv1{2,3},MOVE_base); % convolution
% task regressor: reward delivery
dmFull_task{1,1} = 'task'; % regressor type
dmFull_task{1,2} = 'rwd';  % regressor name
dmFull_task{1,3}(:,1) = histcounts([jkvt(~rwI).trEnd]+1000,time_bin)'; % not-rewarded
dmFull_task{1,3}(:,2) = histcounts([jkvt(rwI).rewardT],time_bin)'; % rewarded
dmFull_task{1,4} = basis.conv(dmFull_task{1,3},TASK_base);

% reach phase hand velocity
[~,k] = func.evtKinematics( jkvt, sessionDur );
hvcmS = k.handVel; % hand velocity in cm/S in ms
hvcmSb = hvcmS(1,time_bin(2:end)); %intm(hvcmS,length(time_bin)-1); % DO NOT USE 'intm' here!

hvcmSbReach = smooth2a(max(0,hvcmSb),0,3); % reach phase (away from the initial position)
hvcmSbReachBound = max(hvcmSbReach)*.8; % just use zero-to-max range
[hV_ReachBase, hV_ReachBins, hV_ReachFunc] = basis.linear_cos(p.N_HANDVEL, [0 hvcmSbReachBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to
X_reachVel = hV_ReachFunc(hvcmSbReach); % convolution

% pull phase hand velocity
absHvcmSbPull = smooth2a(abs(min(0,hvcmSb)),0,3); % pull phase (towards the initial position)
absHvcmSbPullBound = max(absHvcmSbPull)*.8; % just use zero-to-absMax range
[hV_PullBase, hV_PullBins, hV_PullFunc] = basis.linear_cos(p.N_HANDVEL, [0 absHvcmSbPullBound], 1, false); % use velocity instead of speed, just edit linear_cos.m to
X_pullVel = hV_PullFunc(absHvcmSbPull); % convolution

dmFull_handVel{1,1} = 'continuous'; % regressor type
dmFull_handVel{1,2} = 'reachVel';   % regressor name
dmFull_handVel{1,3} = hvcmSbReach'; % continuous variable
dmFull_handVel{1,4} = X_reachVel;

dmFull_handVel{2,1} = 'continuous'; % regressor type
dmFull_handVel{2,2} = 'pullVel';    % regressor name
dmFull_handVel{2,3} = absHvcmSbPull'; % continuous variable
dmFull_handVel{2,4} = X_pullVel;

% coarse binning of reach and pull speed
ratio_time = floor(1 / (p.DT/1000)); % ratio of 1000ms to p.DT
ratio_time_100ms = floor(1 / (p.DT/100)); % ratio of 100ms to p.DT
coarse_100msBin = @(x) mean(reshape(x(1:floor(n_bin/ratio_time_100ms)*ratio_time_100ms), ratio_time_100ms, []))';   % mean with 100-ms bin
coarseMax_100msBin = @(x) max(reshape(x(1:floor(n_bin/ratio_time_100ms)*ratio_time_100ms), ratio_time_100ms, []))'; % max with 100-ms bin

reachVel_100ms = coarseMax_100msBin(hvcmSbReach);
reachVel_100msEdge = round(linspace(0,ceil(range(reachVel_100ms)),11)); % to get 10 bins %round(linspace(0,ceil(range(reachVel_100ms)*.8),11)); % to get 10 bins
pullVel_100ms = coarseMax_100msBin(absHvcmSbPull);
pullVel_100msEdge = round(linspace(0,ceil(range(pullVel_100ms)),11)); %round(linspace(0,ceil(range(pullVel_100ms)*.8),11));

% merge and normalize
[X_full, prm_full] = func.normalize_add_constant([dmFull_mv1(:,4)', dmFull_task(1,4)', dmFull_handVel(:,4)']);
tV.X_full = X_full;
tV.prm_full = prm_full;

% organize trials for 5-fold cross validation (type-by-fold)
%allTrials = cell(length(jkvt),1); %allTrials(:,1) = deal({'all'});
[train_mv, test_mv, type_mv] = crossValTrainTestSetsEachTrialType(posTqC, 5, pStartI);
cumsum_pStartI = cumsum(pStartI);
train_mv_pStartI = cellfun(@(a) cumsum_pStartI(a), train_mv, 'un', 0);
test_mv_pStartI  = cellfun(@(a) cumsum_pStartI(a), test_mv, 'un', 0);

rwC = cell(length(rwI),1);
rwC(rwI,1) = deal({'rwd'});
rwC(~rwI,1) = deal({'no'});
[train_rwd,test_rwd,type_rwd] = crossValTrainTestSetsEachTrialType(rwC, 5);

%% 3.1. Design matrix - discrete variables (movement and task regressors)
for f = 1:5 % 5-fold cross-validation for decoding    
    % train trials movement-related variables
    trainLe = false(size(jkvt,2),1); trainLe([train_mv{1,f}; train_mv{2,f}]) = true;
    trainRi = false(size(jkvt,2),1); trainRi([train_mv{3,f}; train_mv{4,f}]) = true;
    trainLo = false(size(jkvt,2),1); trainLo([train_mv{1,f}; train_mv{3,f}]) = true;
    trainHi = false(size(jkvt,2),1); trainHi([train_mv{2,f}; train_mv{4,f}]) = true;
    
    % train trials reward
    trainNoRwd = false(size(jkvt,2),1); trainNoRwd(train_rwd{1,f}) = true;
    trainRwd   = false(size(jkvt,2),1); trainRwd(train_rwd{2,f}) = true;
    
    allTrainThisFold = cell2mat(train_mv_pStartI(:,f));
    allTestThisFold = cell2mat(test_mv_pStartI(:,f));
    allTestThisFold_Rwd = cell2mat(test_rwd(:,f));
    
    % Discrete regressors position
    dm_mv1{1,1} = 'move';
    dm_mv1{1,2} = 'position';
    dm_mv1{1,3}(:,1) = histcounts([jkvt(trainLe).rStartToPull]-2000,time_bin)'; % left success trials
    dm_mv1{1,3}(:,2) = histcounts([jkvt(trainRi).rStartToPull]-2000,time_bin)'; % right success trials
    dm_mv1{1,4} = basis.conv(dm_mv1{1,3},MOVE_base); % convolution
    % Discrete regressors torque
    dm_mv1{2,1} = 'move';
    dm_mv1{2,2} = 'torque';
    dm_mv1{2,3}(:,1) = histcounts([jkvt(trainLo).rStartToPull]-2000,time_bin)'; % low torque success trials
    dm_mv1{2,3}(:,2) = histcounts([jkvt(trainHi).rStartToPull]-2000,time_bin)'; % high torque success trials
    dm_mv1{2,4} = basis.conv(dm_mv1{2,3},MOVE_base); % convolution
    % task regressor: reward delivery
    dm_task{1,1} = 'task'; % regressor type
    dm_task{1,2} = 'rwd';  % regressor name
    dm_task{1,3}(:,1) = histcounts([jkvt(trainNoRwd).trEnd]+1000,time_bin)'; % not-rewarded
    dm_task{1,3}(:,2) = histcounts([jkvt(trainRwd).rewardT],time_bin)'; % rewarded
    dm_task{1,4} = basis.conv(dm_task{1,3},TASK_base);
    
    %% 3.2. Design matrix - continuous variables hand velocity variables (reach and pull phase absolute speeds separately
    % merge and normalize
    [X, prm] = func.normalize_add_constant([dm_mv1(:,4)', dm_task(1,4)', dmFull_handVel(:,4)']);
    tV.X = X;
    tV.prm = prm;
    
    % exclulde segments corresponding to all test sets
    testPtsToExclude = sort([tV.evtRstart(allTestThisFold) tV.evtRwd(allTestThisFold_Rwd)]);
    [~,~,testPtsToExclude_bin] = histcounts(testPtsToExclude,time_bin);
    testDatI = excludeSegments(X,testPtsToExclude_bin,2000./p.DT);
    
    %% 4.1 Firing rates by task variables for initialization of discrete regressors
    %Direction-tuning: 79, 80, 35, 77, 85
    %[i_cell] = getUnitIdSpkTimesCell(spkTimesCell, S, 102); % get the unit index as in the 'spkTimesCell' WR40_081919, unit #102 (#64 in 'spkTimesCell') is an example M1 neuron with target direction tuning
    for i_cell = 1:size(spkTimesCell,2) %1:size(spkTimesCell,2)
        spike_time = clip(spkTimesCell{1,i_cell});
        n_spike = length(spike_time);
        [spike_bin,~,spike_binT] = histcounts(spike_time, time_bin);
        spike_rate = sum(spike_bin) / (sessionDur/(1000)); % spikes/Sec
        spike_bin_conv = basis.normal_filter((spike_bin.*(1000/p.DT))', p.filter_sigma, p.DT); % p.filter_sigma = 100 ms
        
        if spike_rate>=0.5
            % align to the task event (reach start to pull)
            [spike_rStart_count, time_tr_org] = alignToTaskEvent(tV.evtRstart-2000, spike_time, p.DT, p.DT, p.window + 4 * p.filter_sigma);
            spike_rStart = spike_rStart_count'.*(1000./p.DT); % spike rate (Hz)
            in_t = 1 + p.cut:length(time_tr_org) - p.cut;
            time_tr = time_tr_org(in_t);
            
            % binned spike counts aligned to reward time (trEnd+1000)
            [spike_rwdOrNo_count] = alignToTaskEvent(tV.evtRwd, spike_time, p.DT, p.DT, p.window + 4 * p.filter_sigma);
            spike_rwdOrNo = spike_rwdOrNo_count'.*(1000./p.DT); % spike rate (Hz)
            
            spike_leRi = func.group_stat2(spike_rStart, riI(pStartI)+1); % mean PSTH, lelo vs. rest
            spike_loHi = func.group_stat2(spike_rStart, hiI(pStartI)+1); % mean PSTH, lelo vs. rest
            spike_reward = func.group_stat2(spike_rwdOrNo, rwI+1); % mean PSTH, no Reward vs. reward
            
            % normal filter
            spike_leRi_conv = basis.normal_filter(spike_leRi, p.filter_sigma, p.DT);
            spike_loHi_conv = basis.normal_filter(spike_loHi, p.filter_sigma, p.DT);
            spike_reward_conv = basis.normal_filter(spike_reward, p.filter_sigma, p.DT);
            
            % log
            leRi_log = log(spike_leRi_conv(in_t, :)/spike_rate);
            loHi_log = log(spike_loHi_conv(in_t, :)/spike_rate);
            reward_log = log(spike_reward_conv(in_t,:)/spike_rate);
            
            %% 4.2 Firing rates by task variables for initialization of continuous regressors
            spike_100ms = coarse_100msBin(spike_bin) * ratio_time; % spike rate 100-ms bin in Hz (convertion to Hz - multiply the ratio of 1000ms to the original binSize)
            
            % bin max reach velocity
            [rVel_mSpkR_100ms, rVel_sSpkR_100ms, rVel_bin_100ms] = func.group_stat(reachVel_100ms, spike_100ms, reachVel_100msEdge); % get the mean and sem spike rates per binned reach velocity
            rVel_logSpkRate_100ms = log(rVel_mSpkR_100ms / spike_rate); %
            
            % bin max pull velocity
            [pVel_mSpkR_100ms, pVel_sSpkR_100ms, pVel_bin_100ms] = func.group_stat(pullVel_100ms, spike_100ms, pullVel_100msEdge); % get the mean and sem spike rates per binned reach velocity
            pVel_logSpkRate_100ms = log(pVel_mSpkR_100ms / spike_rate);
            
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
            w_c0 = log(spike_rate);
            
            % rebuilding fitted kernel
            leRi0 = rStart_base * w_leRi0;
            loHi0 = rStart_base * w_loHi0;
            reward0 = reward_base * w_reward0;
            rVel0 = rVel_base * w_rVel0; % reach velocity
            pVel0 = pVel_base * w_pVel0; % pull velocity
            
            %% 6. loss function & optimization
            algopts = {'algorithm','trust-region','Gradobj','on','Hessian','on', 'display', 'iter', 'maxiter', 100};
            opts = optimset(algopts{:});
            % initialize weights
            w0s = [w_c0; w_leRi0(:); w_loHi0(:); w_reward0(:); w_rVel0; w_pVel0];
            if CALC_PRM && sum(isnan(w0s))==0 && sum(isinf(w0s))==0
                prm0 = [w_c0; w_leRi0(:); w_loHi0(:); w_reward0(:); w_rVel0; w_pVel0];
            else
                w_c0 = log(spike_rate);
                prm0 = [w_c0; rand(prm.n_var, 1) - 0.5];
            end
            
            prm0_norm = [prm0(1) + prm.mean * prm0(2:end); prm0(2:end) .* prm.std'];
            
            % loss functions (full, drop, single)
            lfunc = @(w) loss.log_poisson_loss(w, X(~testDatI,:), spike_bin(~testDatI)', p.DT/1000); % make sure that p.DT here should be 1/1000, well not necessarily it depends on the binSize
            
            % get cvR2 and dR2 with cross-validation
            for rr = 1:length(prm.index) % note that prm.index{1} contains ones (constant)
                if rr == 1 % fit without dropping
                    % fit using all regressor sets without dropping
                    [prm1_norm, loss1, exitflag, output, grad, hessian] = fminunc(lfunc, prm0_norm, opts); % at 1st iteration, fit with all kernels
                    Xprm1_full = X_full*prm1_norm;  % make sure to use 'X_full' when evaluating R^2 of the test sets
                    expXprm1_full = basis.normal_filter(exp(Xprm1_full), p.filter_sigma, p.DT);
                    % calculate r^2 within PSTHs
                    rez(i_cell).cvR2_psth_train{f,rr} = getR2psth(tV.evtRstart(allTrainThisFold), time_bin, expXprm1_full, spike_bin_conv, 10, 2000);
                    rez(i_cell).cvR2_psth_test{f,rr} = getR2psth(tV.evtRstart(allTestThisFold), time_bin, expXprm1_full, spike_bin_conv, 10, 2000);
                    rez(i_cell).dR2_psth_train{f,rr} = NaN; 
                    rez(i_cell).dR2_psth_test{f,rr} = NaN;
                    
                elseif rr >= 2
                    % fit dropping each regressor set for dR2
                    prmDropI = ~ismember(1:length(prm0),prm.index{rr});  % include all except the current set
                    prm0_norm_drop = prm0_norm(prmDropI);
                    X_drop = X(:,prmDropI);
                    lfunc_drop   = @(w) loss.log_poisson_loss(w, X_drop(~testDatI,:), spike_bin(~testDatI)', p.DT/1000);
                    [prm1_norm_drop, ~, ~, ~, ~, ~] = fminunc(lfunc_drop, prm0_norm_drop, opts); % fit
                    Xprm1_drop_full = X_full(:,prmDropI)*prm1_norm_drop;  % make sure to use 'X_full' when evaluating R^2 of the test sets
                    expXprm1_drop_full = basis.normal_filter(exp(Xprm1_drop_full), p.filter_sigma, p.DT);
                    tempTrainR2 = getR2psth(tV.evtRstart(allTrainThisFold), time_bin, expXprm1_drop_full, spike_bin_conv, 10, 2000);
                    tempTestR2 = getR2psth(tV.evtRstart(allTestThisFold), time_bin, expXprm1_drop_full, spike_bin_conv, 10, 2000);
                    rez(i_cell).dR2_psth_train{f,rr} = rez(i_cell).cvR2_psth_train{f,1}-tempTrainR2; % get dR2 relative to the full model fit
                    rez(i_cell).dR2_psth_test{f,rr} = rez(i_cell).cvR2_psth_test{f,1}-tempTestR2; % get dR2 relative to the full model fit
                    
                    % fit with each regressor set only for cvR2
                    prmSingI = ismember(1:length(prm0),prm.index{rr}); % include only the current set
                    prmSingI(1) = true; % to keep the constant
                    prm0_norm_sing = prm0_norm(prmSingI);
                    X_sing = X(:,prmSingI);
                    lfunc_sing = @(w) loss.log_poisson_loss(w, X_sing(~testDatI,:), spike_bin(~testDatI)', p.DT/1000);
                    [prm1_norm_sing, ~, ~, ~, ~, ~] = fminunc(lfunc_sing, prm0_norm_sing, opts); % fit
                    Xprm1_sing_full = X_full(:,prmSingI)*prm1_norm_sing;  % make sure to use 'X_full' when evaluating R^2 of the test sets
                    expXprm1_sing_full = basis.normal_filter(exp(Xprm1_sing_full), p.filter_sigma, p.DT);
                    rez(i_cell).cvR2_psth_train{f,rr} = getR2psth(tV.evtRstart(allTrainThisFold), time_bin, expXprm1_sing_full, spike_bin_conv, 10, 2000);
                    rez(i_cell).cvR2_psth_test{f,rr} = getR2psth(tV.evtRstart(allTestThisFold), time_bin, expXprm1_sing_full, spike_bin_conv, 10, 2000);                  
                end
            end
            fprintf('processed cell # %d for fold # %d\n', i_cell, f) % report unit progression
        end
    end
end


%save(fullfile(filePath,strcat('glm_cvR2_dR2_',saveName)),'rez','tV','p')
%save(fullfile(filePath, strcat('glm_cvR2_dR2_',saveName),'rez','tV','p')
    

    %% individual unit psth aligned to a task event
    % % sort trials
    % tq = round([jkvt.pull_torque]./10); % torque
    % tqT = tq(S.trI)'; tqT(:,2) = 1:length(S.trI); tqT(:,3) = S.trI;  sortByTq = sortrows(tqT,1);
    % ps = [jkvt.reachP1]; % position
    % psT = ps(S.trI)'; psT(:,2) = 1:length(S.trI); psT(:,3) = S.trI;  sortByPs = sortrows(psT,1);
    % tqpsT(:,1) = tqT(:,1)+psT(:,1); tqpsT(:,2) = 1:length(S.trI); tqpsT(:,3) = S.trI; sortByTqPs = sortrows(tqpsT,1); % torque position combo
    %
    % %load('binSpkCountSTRCTXWR40_082019.mat', 'rStartToPull')
    % %S = rStartToPull
    % cellI = 23; % 102 79, 80, 35, 77, 85
    % thisUnitSpkTimes = S.SpkTimes{cellI};
    % individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTq, cellI, [3e3 2e3], [2e3 1.9e3]);
    % individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByPs, cellI, [3e3 2e3], [2e3 1.9e3]);
    % individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTqPs, cellI, [3e3 2e3], [2e3 1.9e3]);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [binSpkMatOut,binEdgesOut] = alignToTaskEvent(evt, spikeTime, binSize, stepSize, oneWindow)
    binEdgesOut = -oneWindow:binSize:oneWindow;
    bin1msEdge = cellfun(@(a) max(a-oneWindow,0):a+oneWindow+1, num2cell(evt),'un', 0);
    bin1msCount = cellfun(@(a) histcounts(spikeTime,a), bin1msEdge, 'un',0);
    binSpkMatOut = bin1msSpkCountMat(cell2mat(bin1msCount'),binSize,stepSize);
    end
    
    function [unitIdSTC] = getUnitIdSpkTimesCell(spkTimesCellIn, sIn, unitIdS)
    %unit Ids might not match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
    % this function identifies the unit from 'S' in 'spkTimesCell' and outputs
    % it as 'unitidSTC'
    
    geomI = cell2mat(cellfun(@(a) isequal(a,sIn.geometry{unitIdS}), spkTimesCellIn(4,:), 'un', 0));
    wfI = cell2mat(cellfun(@(a) isequal(a,sIn.meanWF{unitIdS}), spkTimesCellIn(6,:), 'un', 0));
    
    unitIdSTC = find(geomI & wfI);
    
    end
    
    function [unitIdS] = getUnitIdSfromSpkTimesCell(spkTimesCellIn, sIn, unitIdC)
    %unit Ids might not match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
    % this function identifies the unit from 'S' in 'spkTimesCell' and outputs
    % it as 'unitidSTC'
    
    geomI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{4,unitIdC}), sIn.geometry, 'un', 0));
    wfI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{6,unitIdC}), sIn.meanWF, 'un', 0));
    
    unitIdS = find(geomI & wfI);
    
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
    
    function [trainC,testC,type] = crossValTrainTestSetsEachTrialType(typeC, fold, varargin)
    %This function takes trial type information (e.g. four types (two position-by-two torques)
    % and generates N-fold cross-validation train and test sets. An additional
    % logic can be input for further selection of trials as a varargin.
    % use e.g.1: [trainC,testC] = crossValTrainTestSetsEachTrialType(posTqC, 5, pStartI)
    % use e.g.2:
    
    useAllTrainTrs = true;
    
    if nargin == 2
        addLogic = true(size(typeC,1),1);
    elseif nargin == 3
        addLogic = varargin{1};
    end
    
    assert(length(addLogic)==length(typeC))
    
    type = unique(typeC);
    for ty = 1:length(type)
        typeS{1,ty} = find(cell2mat(cellfun(@(a) strcmpi(a,type{ty}) , typeC, 'un', 0)) & addLogic);
        typeS{1,ty} = randsample(typeS{1,ty}, length(typeS{1,ty})); % randomize the order
    end
    
    trialN = min(cellfun(@length,typeS));
    testN = floor(trialN*(1/fold));
    trainN = floor(trialN*((fold-1)/fold));
    
    testC  = cell(length(type),fold); % test sets per type and fold
    trainC = cell(length(type),fold); % train sets per type and fold
    
    for ff = 1:fold
        tempTest = (ff-1)*testN+1:(ff-1)*testN+testN; % test trials of this fold (common across trial types)
        for tp = 1:length(type)
            if useAllTrainTrs
                tempTrain = find(~ismember(1:length(typeS{tp}),tempTest)); % sample train trials of this trial type of this fold
            else
                tempTrain = randsample(find(~ismember(1:length(typeS{tp}),tempTest)),trainN); % sample train trials of this trial type of this fold
            end
            testC{tp,ff}  = sort(typeS{tp}(tempTest));
            trainC{tp,ff} = sort(typeS{tp}(tempTrain));
        end
    end
    end
    
    function [testDatC] = getTestDatCell(dataSet,indices)
    for ii = 1:length(indices)
        testDatC{ii} = dataSet(:,indices(ii));
    end
    end
    
    % compute the probability of being a left-target trial
    function [prob1, post1, post2] = posteriorProb2(model1_lambda, model2_lambda, testDat)
    post1 = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(testDat), model1_lambda, 'un', 0));
    post2 = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(testDat), model2_lambda, 'un', 0));
    prob1 = post1./(post1+post2);
    end
    
    function [excludeLogic] = excludeSegments(dat,pointsToExclude,oneWindow)
    % this function takes a whole segment and excludes
    %subsegment-windows specified by the points and window.
    % dat = X;
    % pointsToExclude = testPtsToExclude_bin;
    % oneWindow = 2000/p.DT;
    
    if size(dat,1)<size(dat,2)
        dat = dat';
    end
    
    windowsToExclude = arrayfun(@(a) a-oneWindow:a+oneWindow, pointsToExclude, 'un', 0);
    excludeLogicC = cellfun(@(a) ismember(1:size(dat,1),a), windowsToExclude, 'un', 0)';
    excludeLogic = sum(cell2mat(excludeLogicC))>0;
    
    if size(excludeLogic,1)<size(excludeLogic,2)
        excludeLogic = excludeLogic';
    end
    
    end
    
    function [r2_psth] = getR2psth(evtInMs, timeBin, rateFit, rateDat, binSize, oneWin)
    
    rr2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
    model_psth = func.alignGlmOutToTaskEvent(evtInMs, timeBin, rateFit, binSize, oneWin); % fitted rate aligned to reachStart
    spike_psth = func.alignGlmOutToTaskEvent(evtInMs, timeBin, rateDat, binSize, oneWin); % fitted rate aligned to reachStart
    r2_psth = max(rr2(reshape(spike_psth,[],1),reshape(model_psth,[],1)),0);
    
    end



