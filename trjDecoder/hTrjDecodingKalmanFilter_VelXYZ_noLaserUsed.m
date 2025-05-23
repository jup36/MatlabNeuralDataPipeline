function hTrjDecodingKalmanFilter_VelXYZ_noLaserUsed(filePath, saveName, plotlogic)
%This decodes kinematics of mouse 3-d hand movement trajectories (X,Y,Z)
% using cross-validated (leave-a-trial-out) Kalman filter decoding.

%plot logic = 0;
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)

kfDir = dir('preprocessKFdecodeHTrjCtxStr_hTrjF*');
load(fullfile(kfDir.folder,kfDir.name),'s')
%fileName = 'preprocessKFdecodeHTrjCtxStr_WR40_081919.mat';
%load(fullfile(filePath,fileName),'s')

figSaveDir = fullfile(filePath,'Figure','KalmanFilter_decoding');
if ~(isfolder(fullfile(filePath,'Figure','KalmanFilter_decoding')))
    mkdir(fullfile(filePath,'Figure','KalmanFilter_decoding'))
end

resample = 100;
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtx, 'un', 0));
if isfield(s.dat,'laserIdx')
    stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, s.dat.laserIdx, 'un', 0)); % stim trials
    trainTrN = min(sum(valTrI & ~stmTrI))-1;
else
    trainTrN = min(sum(valTrI))-1;
end
ctxCnumb = max(unique(cell2mat(cellfun(@(a) size(a,1), s.dat.spkCtx, 'un', 0)))); % # of cortex cells
strCnumb = max(unique(cell2mat(cellfun(@(a) size(a,1), s.dat.spkStr, 'un', 0)))); % # of striatal cells
minCnumb = min(ctxCnumb, strCnumb); % # of cells to be included

%% select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
for r = 1:size(s.dat.state,1)
    for c = 1:size(s.dat.state,2)
        if ~isempty(s.dat.state{r,c})
            s.dat.state{r,c} = s.dat.state{r,c}(4:6,:); % X, Y, Z velocity (row 1:3 for position, row 4:6 for velocity)
        end
    end
end
clearvars r c

%% leave-a-trial-out decoding using Kalman Filter (heavy-lifting part)
for i = 1:resample %resample % repeat resampling trials
    randCtxI = randperm(ctxCnumb);
    randStrI = randperm(strCnumb);
    ctxI = randCtxI(1:minCnumb); % ctx cells for this iteration
    strI = randStrI(1:minCnumb); % str cells for this iteration
    
    for r = 1:size(s.dat.spkCtx,1) % row: trials
        for c = 1:size(s.dat.spkCtx,2) % column: position/torque pairs
            trainI = true(size(valTrI,1),size(valTrI,2));
            if valTrI(r,c)
                trainI(r,c) = false; % to leave one trial out as a test trial
                testI = ~trainI;     % index for the one test trial left out
                
                if isfield(s.dat,'laserIdx')
                    valTrainI = trainI & valTrI & ~stmTrI; % valid train trial index (exclude stim trials)
                else
                    valTrainI = trainI & valTrI; % valid train trial index (when there's no stim trials)
                end
                
                % get current train trials by resampling (to include the same # of trials for each trial type)
                currTrainTrs = zeros(trainTrN,size(trainI,2));
                currTrainState = []; % current training data states (hand position: row 1-3, velocity: row 4-6)
                currTrainCtx = [];
                currTrainStr = [];
                
                for cc = 1:size(valTrainI,2) % sample trials from each trial type
                    tempValTr = find(valTrainI(:,cc));
                    tempValTrRand = tempValTr(randperm(length(tempValTr)));
                    currTrainTrs(:,cc) = tempValTrRand(1:trainTrN); % take the set # of randomized trials from each trial type (column)
                    currTrainState = [currTrainState, s.dat.state(currTrainTrs(:,cc),cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                    %currTrainState = cellfun(@(a) a(1:3,:), currTrainState, 'un', 0); % choose variables to estimate
                    currTrainCtx = [currTrainCtx, cellfun(@(a) a(ctxI,:), s.dat.spkCtx(currTrainTrs(:,cc),cc), 'un', 0)]; % ctx spike mat with cells resampled to match # of cells
                    currTrainStr = [currTrainStr, cellfun(@(a) a(strI,:), s.dat.spkStr(currTrainTrs(:,cc),cc), 'un', 0)]; % str spike mat with cells resampled to match # of cells
                end
                clearvars cc
                [NTrain,NTrType] = size(currTrainState);
                
                %% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                interStateM = 0;    % inter-state matrix
                intraStateM = 0;    % intra-state matrix
                obsStateCtxM = 0;   % observed Ctx state matrix
                obsStateStrM = 0;   % observed Str state matrix
                obsIntraStateM = 0; % observed intra-state matrix
                count = 0;
                countPool = 0;
                for t = 1:size(currTrainState,1) % # of trial
                    for tt = 1:size(currTrainState,2) % # of trial types
                        % for parameter A, state model describes how state evolves over time
                        stateZ1 = currTrainState{t,tt}(:,2:end); % Zt
                        stateZ2 = currTrainState{t,tt}(:,1:end-1); % Zt-1
                        tmp1 = stateZ1*stateZ2'; % Zt*Zt-1' (6-by-6 matrix)
                        interStateM = interStateM+tmp1; % sum Zt*Zt-1'
                        tmp2 = stateZ2*stateZ2'; % Zt-1*Zt-1' (6-by-6 matrix)
                        intraStateM = intraStateM+tmp2; % sum Zt-1*Zt-1'
                        % for parameter C, observation model describes how observation relates to the state
                        stateZ = currTrainState{t,tt}; % state (position, velocity)
                        obsCtx = currTrainCtx{t,tt}; % ctx spike data
                        obsStr = currTrainStr{t,tt}; % str spike data
                        tmp3 = obsCtx*stateZ'; % Xt*Zt' (Ctx cell-by-state mat)
                        obsStateCtxM = obsStateCtxM+tmp3; % sum Ctx observation-state matrix
                        tmp4 = obsStr*stateZ'; % Xt*Zt' (Str cell-by-state mat)
                        obsStateStrM = obsStateStrM+tmp4; % sum Str observation-state matrix
                        tmp5 = stateZ*stateZ'; % Zt*Zt'
                        obsIntraStateM = obsIntraStateM+tmp5; % sum Zt*Zt'
                        % for parameter Pi and V
                        count=count+1;
                        poolZStart(:,count) = currTrainState{t,tt}(:,1);
                    end
                end
                clearvars t tt
                A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
                C_ctx=(obsStateCtxM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx spikes, which defines how observation relates to the state
                C_str=(obsStateStrM/count)/(obsIntraStateM/count); % Slope for the observation model for str spikes
                Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
                V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
                clearvars poolZStart
                
                % Fit Q and R (variance)
                sumQM = 0;
                sumRM_ctx = 0;
                sumRM_str = 0;
                count = 0;
                for t = 1:size(currTrainState,1) % # of trial
                    for tt = 1:size(currTrainState,2) % # of trial types
                        % for parameter Q, covariance of the state model
                        stateZ1 = currTrainState{t,tt}(:,2:end); % Zt
                        stateZ2 = currTrainState{t,tt}(:,1:end-1); % Zt-1
                        tmp6 = stateZ1-A*stateZ2; % Zt-A*Zt-1
                        tmp6 = tmp6*tmp6'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
                        count = count + size(stateZ1,2);
                        sumQM = sumQM + tmp6; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
                        % for parameter R of Ctx spikes, covariance of the observation model
                        stateZ = currTrainState{t,tt}; % Zt
                        obsCtx = currTrainCtx{t,tt}; % ctx spike data
                        tmp7 = obsCtx-C_ctx*stateZ; % Xt_ctx-C_ctx*Zt
                        tmp7 = tmp7*tmp7'; % (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                        sumRM_ctx = sumRM_ctx+tmp7; % Sum (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                        % for parameter R of Str spikes, covariance of the observation model
                        obsStr = currTrainStr{t,tt}; % str spike data
                        tmp8 = obsStr-C_str*stateZ; % Xt_str-C_str*Zt
                        tmp8 = tmp8*tmp8'; % (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                        sumRM_str = sumRM_str+tmp8; % Sum (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    end
                end
                clearvars t tt
                Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
                R_ctx = sumRM_ctx/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                R_str = sumRM_str/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                
                %% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                curTrLength = size(s.dat.state{testI},2); % position, velocity in 20ms bins
                % Initialization
                mu_ctx=Pi;   % sample mean
                sigma_ctx=V; % sample covariance
                mu_str=Pi;   % sample mean
                sigma_str=V; % sample covariance
                clearvars curEstStateMean*;
                clearvars curEstStateCov*;
                for b = 1:curTrLength % # of timebins
                    % One-step prediction by ctx and str data separately
                    mu_ctx = A*mu_ctx; % A is the coefficient matrix for the state model that maps the previous states to current states
                    sigma_ctx = A*sigma_ctx*A'+Q; % Q is the covariance of the state model
                    mu_str = A*mu_str;
                    sigma_str = A*sigma_str*A'+Q;
                    % compute the Kalman gain (needs to separately computed for ctx and str)
                    K_ctx = sigma_ctx*C_ctx'/(C_ctx*sigma_ctx*C_ctx'+R_ctx); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model
                    K_str = sigma_str*C_str'/(C_str*sigma_str*C_str'+R_str); % *inv(C_str*sigma*C_str'+R_str); % Kalman gain for str observation model
                    % update the state by ctx observation
                    curObsX_ctx = s.dat.spkCtx{testI}(ctxI,b); % take the current ctx spikes bin-by-bin
                    mu_ctx = mu_ctx + K_ctx*(curObsX_ctx-C_ctx*mu_ctx);
                    sigma_ctx = sigma_ctx - K_ctx*C_ctx*sigma_ctx;
                    curEstStateMean_ctx(:,b)=mu_ctx;
                    curEstStateCov_ctx(:,:,b)=sigma_ctx;
                    
                    % update the state by str observation
                    curObsX_str = s.dat.spkStr{testI}(strI,b); % take the current str spikes bin-by-bin
                    mu_str = mu_str + K_str*(curObsX_str-C_str*mu_str);
                    sigma_str = sigma_str - K_str*C_str*sigma_str;
                    curEstStateMean_str(:,b)=mu_str;
                    curEstStateCov_str(:,:,b)=sigma_str;
                end
                s.dat.estStateCtxMean{r,c,i} = curEstStateMean_ctx;
                s.dat.estStateCtxCov{r,c,i} = curEstStateCov_ctx;
                
                s.dat.estStateStrMean{r,c,i} = curEstStateMean_str;
                s.dat.estStateStrCov{r,c,i} = curEstStateCov_str;
                
            end
        end
    end
    fprintf('finished iteration# %d\n', i)
end
clearvars r c

%% get the # of kinematic variables
valCellI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.state, 'un', 0));
sizeCell = cellfun(@size, s.dat.state, 'un', 0);
valSizeC = cell2mat(sizeCell(valCellI));
nKv = unique(valSizeC(:,1)); % the # of kinematic variables

%% take average across estimated trajectories of resampling, interpolate to match trajectory lengths
for rr = 1:size(s.dat.estStateCtxMean,1) % trials (row)
    for cc = 1:size(s.dat.estStateCtxMean,2) % trial-types (column)
        p1 = find(s.dat.pullIdx{rr,cc}==1,1,'first'); % pull start
        p2 = find(s.dat.pullIdx{rr,cc}==1,1,'last'); % pull stop
        
        % cortex estimated state
        s.dat.stateCtx{rr,cc} = nanmean(cell2mat(s.dat.estStateCtxMean(rr,cc,:)),3); % average across ctx resampled trials
        tmpCtx = s.dat.stateCtx{rr,cc};
        stateCtx{rr,cc} = nan(nKv,100); % original estimated trajectory without interpolation
        if size(tmpCtx,2)>=5
            stateCtx{rr,cc}(:,1:min(size(tmpCtx,2),100)) = tmpCtx(:,1:min(size(tmpCtx,2),100));
            stateCtxInt{rr,cc} = intm(tmpCtx,100); % original estimated trajectory with interpolation
            
            if ~isempty(p1) && ~ isempty(p2) && p2-p1>=1
                prePullCtx{rr,cc} = intm(tmpCtx(:,1:p1),50); % pre-pull trajectory with interpolation
                pullCtx{rr,cc} = intm(tmpCtx(:,p1:p2),50);   % during-pull trajectory with interpolation
                %pstPullCtx{rr,cc} = intm(tmpCtx(:,p2:end),50); % post-pull trajectory with interpolation
            end
        end
        
        % striatum estimated state
        s.dat.stateStr{rr,cc} = nanmean(cell2mat(s.dat.estStateStrMean(rr,cc,:)),3); % average across str resampled trials
        tmpStr = s.dat.stateStr{rr,cc};
        stateStr{rr,cc} = nan(nKv,100); % original estimated trajectory without interpolation
        if size(tmpStr,2)>=5
            stateStr{rr,cc}(:,1:min(size(tmpStr,2),100)) = tmpStr(:,1:min(size(tmpStr,2),100));
            stateStrInt{rr,cc} = intm(tmpStr,100); % original estimated trajectory with interpolation
            
            if ~isempty(p1) && ~ isempty(p2) && p2-p1>=1
                prePullStr{rr,cc} = intm(tmpStr(:,1:p1),50); % pre-pull trajectory with interpolation
                pullStr{rr,cc} = intm(tmpStr(:,p1:p2),50);   % during-pull trajectory with interpolation
                %pstPullStr{rr,cc} = intm(tmpStr(:,p2:end),50); % post-pull trajectory with interpolation
            end
        end
        
        % actual state
        tmpAct = s.dat.state{rr,cc};
        state{rr,cc} = nan(nKv,100);
        
        if size(tmpAct,2)>=5
            state{rr,cc}(:,1:min(size(tmpAct,2),100)) = tmpAct(:,1:min(size(tmpAct,2),100));
            stateInt{rr,cc} = intm(tmpAct,100);
            
            if ~isempty(p1) && ~ isempty(p2) && p2-p1>=1
                prePullAct{rr,cc} = intm(tmpAct(:,1:p1),50); % pre-pull trajectory with interpolation
                pullAct{rr,cc} = intm(tmpAct(:,p1:p2),50);   % during-pull trajectory with interpolation
                %pstPullAct{rr,cc} = intm(tmpStr(:,p2:end),50); % post-pull trajectory with interpolation
            end
        end
        
        % distance between the actual trajectory and the mean estimated cortex and striatum trajectory
        if ~isempty(s.dat.state{rr,cc})
            s.dat.dTrjCtx{rr,cc} = sqrt((s.dat.state{rr,cc}(1:3,:)-s.dat.stateCtx{rr,cc}).^2);
            s.dat.dTrjStr{rr,cc} = sqrt((s.dat.state{rr,cc}(1:3,:)-s.dat.stateStr{rr,cc}).^2);
        else
            s.dat.dTrjCtx{rr,cc} = sqrt((s.dat.state{rr,cc}-s.dat.stateCtx{rr,cc}).^2);
            s.dat.dTrjStr{rr,cc} = sqrt((s.dat.state{rr,cc}-s.dat.stateStr{rr,cc}).^2);
        end
    end
end
clearvars rr cc

%% Get the mean and sem trajectories for each diraction-torque combinations
valCellI = cell2mat(cellfun(@(a) ~isempty(a), stateCtxInt, 'un', 0));
sizeCell = cellfun(@size, stateCtxInt, 'un', 0);
valSizeC = cell2mat(sizeCell(valCellI));
nTb = unique(valSizeC(:,2)); % the # of time bins

for cc = 1:size(s.dat.estStateCtxMean,2) % trial-types
    %% cortex-estimated interpolated trajectories
    % cortex estimated trajectory with interpolation all valid trials
    tmpKvTimeTrialCtxInt = cell2mat(reshape(stateCtxInt(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMCtxInt{1,cc} = nanmean(tmpKvTimeTrialCtxInt,3);
    s.dat.trSCtxInt{1,cc} = nanstd(tmpKvTimeTrialCtxInt,0,3)./sqrt(size(tmpKvTimeTrialCtxInt,3));
    % cortex estimated trajectory with interpolation all valid trials
    tmpKvTimeTrialCtxPrePull = cell2mat(reshape(prePullCtx(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMCtxPrePull{1,cc} = nanmean(tmpKvTimeTrialCtxPrePull,3);
    s.dat.trSCtxPrePull{1,cc} = nanstd(tmpKvTimeTrialCtxPrePull,0,3)./sqrt(size(tmpKvTimeTrialCtxPrePull,3));
    % cortex estimated trajectory with interpolation all valid trials
    tmpKvTimeTrialCtxPull = cell2mat(reshape(pullCtx(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMCtxPull{1,cc} = nanmean(tmpKvTimeTrialCtxPull,3);
    s.dat.trSCtxPull{1,cc} = nanstd(tmpKvTimeTrialCtxPull,0,3)./sqrt(size(tmpKvTimeTrialCtxPull,3));
    
    %% striatum-estimated interpolated trajectories
    % striatum estimated trajectory with interpolation all valid trials
    tmpKvTimeTrialStrInt = cell2mat(reshape(stateStrInt(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMStrInt{1,cc} = nanmean(tmpKvTimeTrialStrInt,3);
    s.dat.trSStrInt{1,cc} = nanstd(tmpKvTimeTrialStrInt,0,3)./sqrt(size(tmpKvTimeTrialStrInt,3));
    % striatum estimated trajectory with interpolation all valid pre-pull trials
    tmpKvTimeTrialStrPrePull = cell2mat(reshape(prePullStr(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMStrPrePull{1,cc} = nanmean(tmpKvTimeTrialStrPrePull,3);
    s.dat.trSStrPrePull{1,cc} = nanstd(tmpKvTimeTrialStrPrePull,0,3)./sqrt(size(tmpKvTimeTrialStrPrePull,3));
    % striatum estimated trajectory with interpolation all valid during-pull trials
    tmpKvTimeTrialStrPull = cell2mat(reshape(pullStr(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMStrPull{1,cc} = nanmean(tmpKvTimeTrialStrPull,3);
    s.dat.trSStrPull{1,cc} = nanstd(tmpKvTimeTrialStrPull,0,3)./sqrt(size(tmpKvTimeTrialStrPull,3));
    
    %% actual interpolated trajectories
    % actual hand trajectory with interpolation all valid trials
    tmpKvTimeTrialActInt = cell2mat(reshape(stateInt(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMActInt{1,cc} = nanmean(tmpKvTimeTrialActInt,3);
    s.dat.trSActInt{1,cc} = nanstd(tmpKvTimeTrialActInt,0,3)./sqrt(size(tmpKvTimeTrialActInt,3));
    % actual hand trajectory with interpolation all valid pre-pull trials
    tmpKvTimeTrialActPrePull = cell2mat(reshape(prePullAct(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMActPrePull{1,cc} = nanmean(tmpKvTimeTrialActPrePull,3);
    s.dat.trSActPrePull{1,cc} = nanstd(tmpKvTimeTrialActPrePull,0,3)./sqrt(size(tmpKvTimeTrialActPrePull,3));
    % actual hand trajectory with interpolation all valid during-pull trials
    tmpKvTimeTrialActPull = cell2mat(reshape(pullAct(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMActPull{1,cc} = nanmean(tmpKvTimeTrialActPull,3);
    s.dat.trSActPull{1,cc} = nanstd(tmpKvTimeTrialActPull,0,3)./sqrt(size(tmpKvTimeTrialActPull,3));
    
    %% uninterpolated trajectories cortex striatum actual
    % cortex estimated trajectory without interpolation all valid trials
    tmpKvTimeTrialCtx = cell2mat(reshape(stateCtx(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMCtx{1,cc} = nanmean(tmpKvTimeTrialCtx,3);
    s.dat.trSCtx{1,cc} = nanstd(tmpKvTimeTrialCtx,0,3)./sqrt(size(tmpKvTimeTrialCtx,3));
    
    % striatum estimated trajectory without interpolation all valid trials
    tmpKvTimeTrialStr = cell2mat(reshape(stateStr(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMStr{1,cc} = nanmean(tmpKvTimeTrialStr,3);
    s.dat.trSStr{1,cc} = nanstd(tmpKvTimeTrialStr,0,3)./sqrt(size(tmpKvTimeTrialStr,3));
    
    % actual hand trajectory without interpolation all valid trials
    tmpKvTimeTrialAct = cell2mat(reshape(state(valCellI(:,cc),cc),1,1,[]));
    s.dat.trMAct{1,cc} = nanmean(tmpKvTimeTrialAct,3);
    s.dat.trSAct{1,cc} = nanstd(tmpKvTimeTrialAct,0,3)./sqrt(size(tmpKvTimeTrialAct,3));
end
clearvars cc
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrVel_',saveName)),'s') % last saved after training without stim trials 5/27 Wed 9pm
%load(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrVel_',saveName)),'s') % last saved after training without stim trials 5/27 Wed 9pm

%% plot X, Y, Z trajectories separately for high versus low torques
colorMap = [[100 149 237]./255; [50 205 50]./255; [50 50 50]./255]; % colorMap for cortex and striatum
%% X trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory
if plotlogic == 1
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,1}(1,:),s.dat.trSCtxInt{1,1}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,1}(1,:),s.dat.trSStrInt{1,1}(1,:), ...
        1:nTb,s.dat.trMActInt{1,1}(1,:),s.dat.trSActInt{1,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_whole_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,1}(1,:),s.dat.trSCtxPrePull{1,1}(1,:), ...
        1:50,s.dat.trMStrPrePull{1,1}(1,:),s.dat.trSStrPrePull{1,1}(1,:), ...
        1:50,s.dat.trMActPrePull{1,1}(1,:),s.dat.trSActPrePull{1,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_prePull_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,1}(1,:),s.dat.trSCtxPull{1,1}(1,:), ...
        1:50,s.dat.trMStrPull{1,1}(1,:),s.dat.trSStrPull{1,1}(1,:), ...
        1:50,s.dat.trMActPull{1,1}(1,:),s.dat.trSActPull{1,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_pull_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% X trj, low torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,3}(1,:),s.dat.trSCtxInt{1,3}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,3}(1,:),s.dat.trSStrInt{1,3}(1,:), ...
        1:nTb,s.dat.trMActInt{1,3}(1,:),s.dat.trSActInt{1,3}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_whole_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,3}(1,:),s.dat.trSCtxPrePull{1,3}(1,:), ...
        1:50,s.dat.trMStrPrePull{1,3}(1,:),s.dat.trSStrPrePull{1,3}(1,:), ...
        1:50,s.dat.trMActPrePull{1,3}(1,:),s.dat.trSActPrePull{1,3}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_prePull_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,3}(1,:),s.dat.trSCtxPull{1,3}(1,:), ...
        1:50,s.dat.trMStrPull{1,3}(1,:),s.dat.trSStrPull{1,3}(1,:), ...
        1:50,s.dat.trMActPull{1,3}(1,:),s.dat.trSActPull{1,3}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_pull_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% X trj, high torque, position 1(left), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,2}(1,:),s.dat.trSCtxInt{1,2}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,2}(1,:),s.dat.trSStrInt{1,2}(1,:), ...
        1:nTb,s.dat.trMActInt{1,2}(1,:),s.dat.trSActInt{1,2}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_whole_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,2}(1,:),s.dat.trSCtxPrePull{1,2}(1,:), ...
        1:50,s.dat.trMStrPrePull{1,2}(1,:),s.dat.trSStrPrePull{1,2}(1,:), ...
        1:50,s.dat.trMActPrePull{1,2}(1,:),s.dat.trSActPrePull{1,2}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_prePull_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,2}(1,:),s.dat.trSCtxPull{1,2}(1,:), ...
        1:50,s.dat.trMStrPull{1,2}(1,:),s.dat.trSStrPull{1,2}(1,:), ...
        1:50,s.dat.trMActPull{1,2}(1,:),s.dat.trSActPull{1,2}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_pull_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% X trj, high torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,4}(1,:),s.dat.trSCtxInt{1,4}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,4}(1,:),s.dat.trSStrInt{1,4}(1,:), ...
        1:nTb,s.dat.trMActInt{1,4}(1,:),s.dat.trSActInt{1,4}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_whole_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;

    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,4}(1,:),s.dat.trSCtxPrePull{1,4}(1,:), ...
        1:50,s.dat.trMStrPrePull{1,4}(1,:),s.dat.trSStrPrePull{1,4}(1,:), ...
        1:50,s.dat.trMActPrePull{1,4}(1,:),s.dat.trSActPrePull{1,4}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_prePull_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,4}(1,:),s.dat.trSCtxPull{1,4}(1,:), ...
        1:50,s.dat.trMStrPull{1,4}(1,:),s.dat.trSStrPull{1,4}(1,:), ...
        1:50,s.dat.trMActPull{1,4}(1,:),s.dat.trSActPull{1,4}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('xVel_interp_pull_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Y trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,1}(2,:),s.dat.trSCtxInt{1,1}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,1}(2,:),s.dat.trSStrInt{1,1}(2,:), ...
        1:nTb,s.dat.trMActInt{1,1}(2,:),s.dat.trSActInt{1,1}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_whole_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,1}(2,:),s.dat.trSCtxPrePull{1,1}(2,:), ...
        1:50,s.dat.trMStrPrePull{1,1}(2,:),s.dat.trSStrPrePull{1,1}(2,:), ...
        1:50,s.dat.trMActPrePull{1,1}(2,:),s.dat.trSActPrePull{1,1}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_prePull_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,1}(2,:),s.dat.trSCtxPull{1,1}(2,:), ...
        1:50,s.dat.trMStrPull{1,1}(2,:),s.dat.trSStrPull{1,1}(2,:), ...
        1:50,s.dat.trMActPull{1,1}(2,:),s.dat.trSActPull{1,1}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_pull_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Y trj, low torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,3}(2,:),s.dat.trSCtxInt{1,3}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,3}(2,:),s.dat.trSStrInt{1,3}(2,:), ...
        1:nTb,s.dat.trMActInt{1,3}(2,:),s.dat.trSActInt{1,3}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_whole_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,3}(2,:),s.dat.trSCtxPrePull{1,3}(2,:), ...
        1:50,s.dat.trMStrPrePull{1,3}(2,:),s.dat.trSStrPrePull{1,3}(2,:), ...
        1:50,s.dat.trMActPrePull{1,3}(2,:),s.dat.trSActPrePull{1,3}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_prePull_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,3}(2,:),s.dat.trSCtxPull{1,3}(2,:), ...
        1:50,s.dat.trMStrPull{1,3}(2,:),s.dat.trSStrPull{1,3}(2,:), ...
        1:50,s.dat.trMActPull{1,3}(2,:),s.dat.trSActPull{1,3}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_pull_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Y trj, high torque, position 1(left), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,2}(2,:),s.dat.trSCtxInt{1,2}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,2}(2,:),s.dat.trSStrInt{1,2}(2,:), ...
        1:nTb,s.dat.trMActInt{1,2}(2,:),s.dat.trSActInt{1,2}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_whole_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,2}(2,:),s.dat.trSCtxPrePull{1,2}(2,:), ...
        1:50,s.dat.trMStrPrePull{1,2}(2,:),s.dat.trSStrPrePull{1,2}(2,:), ...
        1:50,s.dat.trMActPrePull{1,2}(2,:),s.dat.trSActPrePull{1,2}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_prePull_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,2}(2,:),s.dat.trSCtxPull{1,2}(2,:), ...
        1:50,s.dat.trMStrPull{1,2}(2,:),s.dat.trSStrPull{1,2}(2,:), ...
        1:50,s.dat.trMActPull{1,2}(2,:),s.dat.trSActPull{1,2}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_pull_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Y trj, high torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,4}(2,:),s.dat.trSCtxInt{1,4}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,4}(2,:),s.dat.trSStrInt{1,4}(2,:), ...
        1:nTb,s.dat.trMActInt{1,4}(2,:),s.dat.trSActInt{1,4}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_whole_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,4}(2,:),s.dat.trSCtxPrePull{1,4}(2,:), ...
        1:50,s.dat.trMStrPrePull{1,4}(2,:),s.dat.trSStrPrePull{1,4}(2,:), ...
        1:50,s.dat.trMActPrePull{1,4}(2,:),s.dat.trSActPrePull{1,4}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_prePull_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,4}(2,:),s.dat.trSCtxPull{1,4}(2,:), ...
        1:50,s.dat.trMStrPull{1,4}(2,:),s.dat.trSStrPull{1,4}(2,:), ...
        1:50,s.dat.trMActPull{1,4}(2,:),s.dat.trSActPull{1,4}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('yVel_interp_pull_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Z trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,1}(3,:),s.dat.trSCtxInt{1,1}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,1}(3,:),s.dat.trSStrInt{1,1}(3,:), ...
        1:nTb,s.dat.trMActInt{1,1}(3,:),s.dat.trSActInt{1,1}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_whole_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,1}(3,:),s.dat.trSCtxPrePull{1,1}(3,:), ...
        1:50,s.dat.trMStrPrePull{1,1}(3,:),s.dat.trSStrPrePull{1,1}(3,:), ...
        1:50,s.dat.trMActPrePull{1,1}(3,:),s.dat.trSActPrePull{1,1}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_prePull_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,1}(3,:),s.dat.trSCtxPull{1,1}(3,:), ...
        1:50,s.dat.trMStrPull{1,1}(3,:),s.dat.trSStrPull{1,1}(3,:), ...
        1:50,s.dat.trMActPull{1,1}(3,:),s.dat.trSActPull{1,1}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_pull_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Z trj, low torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,3}(3,:),s.dat.trSCtxInt{1,3}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,3}(3,:),s.dat.trSStrInt{1,3}(3,:), ...
        1:nTb,s.dat.trMActInt{1,3}(3,:),s.dat.trSActInt{1,3}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_whole_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,3}(3,:),s.dat.trSCtxPrePull{1,3}(3,:), ...
        1:50,s.dat.trMStrPrePull{1,3}(3,:),s.dat.trSStrPrePull{1,3}(3,:), ...
        1:50,s.dat.trMActPrePull{1,3}(3,:),s.dat.trSActPrePull{1,3}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_prePull_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,3}(3,:),s.dat.trSCtxPull{1,3}(3,:), ...
        1:50,s.dat.trMStrPull{1,3}(3,:),s.dat.trSStrPull{1,3}(3,:), ...
        1:50,s.dat.trMActPull{1,3}(3,:),s.dat.trSActPull{1,3}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_pull_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Z trj, high torque, position 1(left), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,2}(3,:),s.dat.trSCtxInt{1,2}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,2}(3,:),s.dat.trSStrInt{1,2}(3,:), ...
        1:nTb,s.dat.trMActInt{1,2}(3,:),s.dat.trSActInt{1,2}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_whole_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,2}(3,:),s.dat.trSCtxPrePull{1,2}(3,:), ...
        1:50,s.dat.trMStrPrePull{1,2}(3,:),s.dat.trSStrPrePull{1,2}(3,:), ...
        1:50,s.dat.trMActPrePull{1,2}(3,:),s.dat.trSActPrePull{1,2}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_prePull_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,2}(3,:),s.dat.trSCtxPull{1,2}(3,:), ...
        1:50,s.dat.trMStrPull{1,2}(3,:),s.dat.trSStrPull{1,2}(3,:), ...
        1:50,s.dat.trMActPull{1,2}(3,:),s.dat.trSActPull{1,2}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_pull_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% Z trj, high torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,4}(3,:),s.dat.trSCtxInt{1,4}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,4}(3,:),s.dat.trSStrInt{1,4}(3,:), ...
        1:nTb,s.dat.trMActInt{1,4}(3,:),s.dat.trSActInt{1,4}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_whole_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PRE-PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPrePull{1,4}(3,:),s.dat.trSCtxPrePull{1,4}(3,:), ...
        1:50,s.dat.trMStrPrePull{1,4}(3,:),s.dat.trSStrPrePull{1,4}(3,:), ...
        1:50,s.dat.trMActPrePull{1,4}(3,:),s.dat.trSActPrePull{1,4}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_prePull_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Z trj, low torque, position 1(left), Cortex, Striatum, Actual PULL trajectory
    figure;
    boundedline(1:50,s.dat.trMCtxPull{1,4}(3,:),s.dat.trSCtxPull{1,4}(3,:), ...
        1:50,s.dat.trMStrPull{1,4}(3,:),s.dat.trSStrPull{1,4}(3,:), ...
        1:50,s.dat.trMActPull{1,4}(3,:),s.dat.trSActPull{1,4}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    print(fullfile(filePath,'Figure','KalmanFilter_decoding',strcat('zVel_interp_pull_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% plot concatenated
    timeX = 0;
    figure; hold on;
    for c = 1:size(s.dat.estStateCtxMean,2)
        for r = 1:20 % just to include the first block only per trial type
            if ~isempty(s.dat.state{r,c})
                ctxTrjX = s.dat.stateCtx{r,c}(1,:); % Ctx X trj (left-right, horizontal hand position)
                strTrjX = s.dat.stateStr{r,c}(1,:); % Str X trj
                actTrjX = s.dat.state{r,c}(1,:); % actual trj
                
                tempX = timeX+2:timeX+length(actTrjX)+1;
                timeX = timeX+length(actTrjX)+1; % update timeX
                
                if s.dat.pos1{r,c}==min([s.dat.pos1{:}]) % left
                    if s.dat.trq{r,c}==min([s.dat.trq{:}]) % low-torque
                        plot(tempX,actTrjX,'k','LineWidth',1);
                        plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',1);
                        plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',1);
                    elseif s.dat.trq{r,c}==max([s.dat.trq{:}]) % high-torque
                        plot(tempX,actTrjX,'k','LineWidth',2);
                        plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',2);
                        plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',2);
                    end
                elseif s.dat.pos1{r,c}==max([s.dat.pos1{:}]) % right
                    if s.dat.trq{r,c}==min([s.dat.trq{:}]) % low-torque
                        plot(tempX,actTrjX,'k','LineWidth',1);
                        plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',1);
                        plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',1);
                    elseif s.dat.trq{r,c}==max([s.dat.trq{:}]) % high-torque
                        plot(tempX,actTrjX,'k','LineWidth',2);
                        plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',2);
                        plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',2);
                    end
                end
            end
        end
    end
    hold off;
    clearvars r c
end % plotlogic

%% compute correlation
stateCC = cell2mat(reshape(s.dat.state,[],1)')'; % concatenated state
stateCClelt = cell2mat(reshape(s.dat.state(:,1),[],1)')'; % actual left low torque
stateCCleht = cell2mat(reshape(s.dat.state(:,2),[],1)')'; % actual left high torque
stateCCrilt = cell2mat(reshape(s.dat.state(:,3),[],1)')'; % actual right low torque
stateCCriht = cell2mat(reshape(s.dat.state(:,4),[],1)')'; % actual right high torque

stateCtxCC = cell2mat(reshape(s.dat.stateCtx,[],1)')'; % concatenated cortex estimated state
stateCtxCClelt = cell2mat(reshape(s.dat.stateCtx(:,1),[],1)')'; % cortex left low torque
stateCtxCCleht = cell2mat(reshape(s.dat.stateCtx(:,2),[],1)')'; % cortex left high torque
stateCtxCCrilt = cell2mat(reshape(s.dat.stateCtx(:,3),[],1)')'; % cortex right low torque
stateCtxCCriht = cell2mat(reshape(s.dat.stateCtx(:,4),[],1)')'; % cortex right high torque

stateStrCC = cell2mat(reshape(s.dat.stateStr,[],1)')'; % concatenated striatum estimated state
stateStrCClelt = cell2mat(reshape(s.dat.stateStr(:,1),[],1)')'; % striatum left low torque
stateStrCCleht = cell2mat(reshape(s.dat.stateStr(:,2),[],1)')'; % striatum left high torque
stateStrCCrilt = cell2mat(reshape(s.dat.stateStr(:,3),[],1)')'; % striatum right low torque
stateStrCCriht = cell2mat(reshape(s.dat.stateStr(:,4),[],1)')'; % striatum right high torque

stateCC_sm = smooth2a(stateCC,4,0);
stateCtxCC_sm = smooth2a(stateCtxCC,4,0);
stateStrCC_sm = smooth2a(stateStrCC,4,0);

% cortex all trials
[corrRez.rCtx,corrRez.pCtx] = corr(stateCC, stateCtxCC, 'Rows','complete');
[corrRez.rCtx_sm,corrRez.pCtx_sm] = corr(stateCC_sm, stateCtxCC_sm, 'Rows','complete');
% cortex high torque
[corrRez.rCtxHtq,corrRez.pCtxHtq] = corr([stateCCleht; stateCCriht],[stateCtxCCleht; stateCtxCCriht],'Rows','complete');
[corrRez.rCtxHtq_sm,corrRez.pCtxHtq_sm] = corr(smooth2a([stateCCleht; stateCCriht],4,0),smooth2a([stateCtxCCleht; stateCtxCCriht],4,0),'Rows','complete');
% cortex low torque
[corrRez.rCtxLtq,corrRez.pCtxLtq] = corr([stateCClelt; stateCCrilt],[stateCtxCClelt; stateCtxCCrilt],'Rows','complete');
[corrRez.rCtxLtq_sm,corrRez.pCtxLtq_sm] = corr(smooth2a([stateCClelt; stateCCrilt],4,0),smooth2a([stateCtxCClelt; stateCtxCCrilt],4,0),'Rows','complete');
% cortex left target
[corrRez.rCtxLe,corrRez.pCtxLe] = corr([stateCClelt; stateCCleht],[stateCtxCClelt; stateCtxCCleht],'Rows','complete');
[corrRez.rCtxLe_sm,corrRez.pCtxLe_sm] = corr(smooth2a([stateCClelt; stateCCleht],4,0),smooth2a([stateCtxCClelt; stateCtxCCleht],4,0),'Rows','complete');
% cortex right target
[corrRez.rCtxRi,corrRez.pCtxRi] = corr([stateCCrilt; stateCCriht],[stateCtxCCrilt; stateCtxCCriht],'Rows','complete');
[corrRez.rCtxRi_sm,corrRez.pCtxRi_sm] = corr(smooth2a([stateCCrilt; stateCCriht],4,0),smooth2a([stateCtxCCrilt; stateCtxCCriht],4,0),'Rows','complete');

% striatum all trials
[corrRez.rStr,corrRez.pStr] = corr(stateCC, stateStrCC, 'Rows','complete');
[corrRez.rStr_sm,corrRez.pStr_sm] = corr(stateCC_sm, stateStrCC_sm, 'Rows','complete');
% striatum high torque
[corrRez.rStrHtq,corrRez.pStrHtq] = corr([stateCCleht; stateCCriht],[stateStrCCleht; stateStrCCriht],'Rows','complete');
[corrRez.rStrHtq_sm,corrRez.pStrHtq_sm] = corr(smooth2a([stateCCleht; stateCCriht],4,0),smooth2a([stateStrCCleht; stateStrCCriht],4,0),'Rows','complete');
% striatum low torque
[corrRez.rStrLtq,corrRez.pStrLtq] = corr([stateCClelt; stateCCrilt],[stateStrCClelt; stateStrCCrilt],'Rows','complete');
[corrRez.rStrLtq_sm,corrRez.pStrLtq_sm] = corr(smooth2a([stateCClelt; stateCCrilt],4,0),smooth2a([stateStrCClelt; stateStrCCrilt],4,0),'Rows','complete');
% striatum left target
[corrRez.rStrLe,corrRez.pStrLe] = corr([stateCClelt; stateCCleht],[stateStrCClelt; stateStrCCleht],'Rows','complete');
[corrRez.rStrLe_sm,corrRez.pStrLe_sm] = corr(smooth2a([stateCClelt; stateCCleht],4,0),smooth2a([stateStrCClelt; stateStrCCleht],4,0),'Rows','complete');
% striatum right target
[corrRez.rStrRi,corrRez.pStrRi] = corr([stateCCrilt; stateCCriht],[stateStrCCrilt; stateStrCCriht],'Rows','complete');
[corrRez.rStrRi_sm,corrRez.pStrRi_sm] = corr(smooth2a([stateCCrilt; stateCCriht],4,0),smooth2a([stateStrCCrilt; stateStrCCriht],4,0),'Rows','complete');

% %% compute correlation
% % correlation for all time points
% stateCC = cell2mat(reshape(s.dat.state,[],1)')'; % concatenated state
% stateCtxCC = cell2mat(reshape(s.dat.stateCtx,[],1)')'; % concatenated cortex estimated state
% stateStrCC = cell2mat(reshape(s.dat.stateStr,[],1)')'; % concatenated striatum estimated state
% 
% stateCC_sm = smooth2a(stateCC,4,0);
% stateCtxCC_sm = smooth2a(stateCtxCC,4,0);
% stateStrCC_sm = smooth2a(stateStrCC,4,0);
% 
% [corrRez.rCtx,corrRez.pCtx] = corr(stateCC, stateCtxCC, 'Rows','complete');
% [corrRez.rStr,corrRez.pStr] = corr(stateCC, stateStrCC, 'Rows','complete');
% 
% [corrRez.rCtx_sm,corrRez.pCtx_sm] = corr(stateCC_sm, stateCtxCC_sm, 'Rows','complete');
% [corrRez.rStr_sm,corrRez.pStr_sm] = corr(stateCC_sm, stateStrCC_sm, 'Rows','complete');

%% correlation for reach and pull phases separately
pull1C = cellfun(@(a) find(a,1,'first'), s.dat.pullIdx, 'un',0); % pull start points for each trajectory
pull2C = cellfun(@(a) find(a,1,'last'), s.dat.pullIdx, 'un',0); % pull end points for each trajectory

% get reach and pull phase trajectories
for c = 1:size(s.dat.state,2) 
    for r = 1:size(s.dat.state,1)
        if ~isempty(s.dat.state{r,c}) && ~isempty(s.dat.stateCtx{r,c}) && ~isempty(s.dat.stateStr{r,c}) && ~isempty(pull1C{r,c})
           p1 = pull1C{r,c}; 
           p2 = pull2C{r,c}; 
           % reach phase trajectory
           s.dat.stateR{r,c} = s.dat.state{r,c}(:,1:pull1C{r,c});
           s.dat.stateCtxR{r,c} = s.dat.stateCtx{r,c}(:,1:pull1C{r,c});
           s.dat.stateStrR{r,c} = s.dat.stateStr{r,c}(:,1:pull1C{r,c});
           
           % reach endPoint offset
           s.dat.rEndOffCtx{r,c} = min(abs(s.dat.state{r,c}(:,p1)-s.dat.stateCtx{r,c}(:,p1-2:p1)),[],2);
           s.dat.rEndOffStr{r,c} = min(abs(s.dat.state{r,c}(:,p1)-s.dat.stateStr{r,c}(:,p1-2:p1)),[],2);
           
           % pull phase trajectory 
           s.dat.stateP{r,c} = s.dat.state{r,c}(:,pull1C{r,c}:end);
           s.dat.stateCtxP{r,c} = s.dat.stateCtx{r,c}(:,pull1C{r,c}:end);
           s.dat.stateStrP{r,c} = s.dat.stateStr{r,c}(:,pull1C{r,c}:end);
        
           % pull endPoint offset 
           s.dat.pEndOffCtx{r,c} = min(abs(s.dat.state{r,c}(:,p2)-s.dat.stateCtx{r,c}(:,p2-2:p2+1)),[],2);
           s.dat.pEndOffStr{r,c} = min(abs(s.dat.state{r,c}(:,p2)-s.dat.stateStr{r,c}(:,p2-2:p2+1)),[],2);     
        end
    end
end
clearvars r c 

[trjOffset.mREndOffCtx,~,trjOffset.sREndOffCtx] = meanstdsem(cell2mat(reshape(s.dat.rEndOffCtx,[],1)')'); 
[trjOffset.mREndOffStr,~,trjOffset.sREndOffStr] = meanstdsem(cell2mat(reshape(s.dat.rEndOffStr,[],1)')'); % concatenated state

[trjOffset.mPEndOffCtx,~,trjOffset.sPEndOffCtx] = meanstdsem(cell2mat(reshape(s.dat.pEndOffCtx,[],1)')'); 
[trjOffset.mPEndOffStr,~,trjOffset.sPEndOffStr] = meanstdsem(cell2mat(reshape(s.dat.pEndOffStr,[],1)')'); % concatenated state

% reach phase correlation
stateCcR = cell2mat(reshape(s.dat.stateR,[],1)')'; % concatenated state
stateCtxCcR = cell2mat(reshape(s.dat.stateCtxR,[],1)')'; % concatenated cortex estimated state
stateStrCcR = cell2mat(reshape(s.dat.stateStrR,[],1)')'; % concatenated striatum estimated state

stateCcR_sm = smooth2a(stateCcR,4,0);
stateCtxCcR_sm = smooth2a(stateCtxCcR,4,0);
stateStrCcR_sm = smooth2a(stateStrCcR,4,0);

[corrRez.rCtxRch,corrRez.pCtxRch] = corr(stateCcR, stateCtxCcR, 'Rows','complete');
[corrRez.rStrRch,corrRez.pStrRch] = corr(stateCcR, stateStrCcR, 'Rows','complete');

[corrRez.rCtxRch_sm,corrRez.pCtxRch_sm] = corr(stateCcR_sm, stateCtxCcR_sm, 'Rows','complete');
[corrRez.rStrRch_sm,corrRez.pStrRch_sm] = corr(stateCcR_sm, stateStrCcR_sm, 'Rows','complete');

% pull phase correlation
stateCcP = cell2mat(reshape(s.dat.stateP,[],1)')'; % concatenated state
stateCtxCcP = cell2mat(reshape(s.dat.stateCtxP,[],1)')'; % concatenated cortex estimated state
stateStrCcP = cell2mat(reshape(s.dat.stateStrP,[],1)')'; % concatenated striatum estimated state

stateCcP_sm = smooth2a(stateCcP,4,0);
stateCtxCcP_sm = smooth2a(stateCtxCcP,4,0);
stateStrCcP_sm = smooth2a(stateStrCcP,4,0);

[corrRez.rCtxPul,corrRez.pCtxPul] = corr(stateCcP, stateCtxCcP, 'Rows','complete');
[corrRez.rStrPul,corrRez.pStrPul] = corr(stateCcP, stateStrCcP, 'Rows','complete');

[corrRez.rCtxPul_sm,corrRez.pCtxPul_sm] = corr(stateCcP_sm, stateCtxCcP_sm, 'Rows','complete');
[corrRez.rStrPul_sm,corrRez.pStrPul_sm] = corr(stateCcP_sm, stateStrCcP_sm, 'Rows','complete');

save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrVel_',saveName)),'corrRez','trjOffset','s','-append')

%% compute residuals errors (1. from the beginning, 2. after pull stop)
% for r = 1:size(s.dat.estStateCtxMean,1)
%     for c = 1:size(s.dat.estStateCtxMean,2)
%         if ~isempty(s.dat.state{r,c})
%             % get pullStop point
%             pStop = find(s.dat.pullIdx{r,c},1,'last');
%             if isempty(pStop)
%                 pStop = size(s.dat.pullIdx{r,c},2);
%             end
%             rsdCtx1{r,c} = s.dat.stateCtx{r,c}(:,1:pStop)-s.dat.state{r,c}(:,1:pStop);
%             rsdCtx1{r,c} = sqrt(rsdCtx1{r,c}.^2); % squared error
%             % interpolation
%             if size(rsdCtx1{r,c},2)>=5
%                 if size(rsdCtx1{r,c},2)>=100
%                     rsdCtx1Int{r,c} = rsdCtx1{r,c}(:,1:100);
%                 elseif size(rsdCtx1{r,c},2)<100
%                     rsdCtx1Int{r,c} = intm(rsdCtx1{r,c},100);
%                 end
%             elseif size(rsdCtx1{r,c},2)<5
%                 rsdCtx1Int{r,c} = nan(size(rsdCtx1{r,c},1),100);
%             end
%             
%             % interpolation
%             rsdCtx2{r,c} = s.dat.stateCtx{r,c}(:,pStop:end)-s.dat.state{r,c}(:,pStop:end);
%             rsdCtx2{r,c} = sqrt(rsdCtx2{r,c}.^2);
%             if size(rsdCtx2{r,c},2)>=5
%                 if size(rsdCtx2{r,c},2)>=100
%                     rsdCtx2Int{r,c} = rsdCtx2{r,c}(:,1:100);
%                 elseif size(rsdCtx2{r,c},2)<100
%                     rsdCtx2Int{r,c} = intm(rsdCtx2{r,c},100);
%                 end
%             elseif size(rsdCtx2{r,c},2)<5
%                 rsdCtx2Int{r,c} = nan(size(rsdCtx2{r,c},1),100);
%             end
%             
%             % interpolation
%             rsdStr1{r,c} = s.dat.stateStr{r,c}(:,1:pStop)-s.dat.state{r,c}(:,1:pStop);
%             rsdStr1{r,c} = sqrt(rsdStr1{r,c}.^2);
%             if size(rsdStr1{r,c},2)>=5
%                 if size(rsdStr1{r,c},2)>=100
%                     rsdStr1Int{r,c} = rsdStr1{r,c}(:,1:100);
%                 elseif size(rsdStr1{r,c},2)<100
%                     rsdStr1Int{r,c} = intm(rsdStr1{r,c},100);
%                 end
%             elseif size(rsdStr1{r,c},2)<5
%                 rsdStr1Int{r,c} = nan(size(rsdStr1{r,c},1),100);
%             end
%             
%             % interpolation
%             rsdStr2{r,c} = s.dat.stateStr{r,c}(:,pStop:end)-s.dat.state{r,c}(:,pStop:end);
%             rsdStr2{r,c} = sqrt(rsdStr2{r,c}.^2);
%             
%             if size(rsdStr2{r,c},2)>=5
%                 if size(rsdStr2{r,c},2)>=100
%                     rsdStr2Int{r,c} = rsdStr2{r,c}(:,1:100);
%                 elseif size(rsdStr2{r,c},2)<100
%                     rsdStr2Int{r,c} = intm(rsdStr2{r,c},100);
%                 end
%             elseif size(rsdStr2{r,c},2)<5
%                 rsdStr2Int{r,c} = nan(size(rsdStr2{r,c},1),100);
%             end
%         end
%     end
% end
% clearvars r c
% 
% msRsdToPullStopCtx = cell(2,size(s.dat.state,2));
% msRsdFrPullStopCtx = cell(2,size(s.dat.state,2));
% msRsdToPullStopStr = cell(2,size(s.dat.state,2));
% msRsdFrPullStopStr = cell(2,size(s.dat.state,2));
% 
% for c = 1:size(rsdCtx1Int,2)
%     
%     currToPullStopCtx = reshape([rsdCtx1Int{:,c}],size(rsdCtx1Int{c},1),size(rsdCtx1Int{c},2),[]);
%     currFrPullStopCtx = reshape([rsdCtx2Int{:,c}],size(rsdCtx2Int{c},1),size(rsdCtx2Int{c},2),[]);
%     
%     currToPullStopStr = reshape([rsdStr1Int{:,1}],size(rsdStr1Int{1},1),size(rsdStr1Int{1},2),[]);
%     currFrPullStopStr = reshape([rsdStr2Int{:,1}],size(rsdStr2Int{1},1),size(rsdStr2Int{1},2),[]);
%     
%     for tt = 1:3 % increament kinematic variables (x,y,z) position
%         [msRsdToPullStopCtx{1,c}(tt,:),~,msRsdToPullStopCtx{2,c}(tt,:)] = meanstdsem(squeeze(currToPullStopCtx(tt,:,:))');
%         [msRsdFrPullStopCtx{1,c}(tt,:),~,msRsdFrPullStopCtx{2,c}(tt,:)] = meanstdsem(squeeze(currFrPullStopCtx(tt,:,:))');
%         
%         [msRsdToPullStopStr{1,c}(tt,:),~,msRsdToPullStopStr{2,c}(tt,:)] = meanstdsem(squeeze(currToPullStopStr(tt,:,:))');
%         [msRsdFrPullStopStr{1,c}(tt,:),~,msRsdFrPullStopStr{2,c}(tt,:)] = meanstdsem(squeeze(currFrPullStopStr(tt,:,:))');
%     end
% end
% clearvars c tt
% 
% %% plot residual
% colorMap = [[100 149 237]./255; [50 205 50]./255];
% figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(1,:), msRsdToPullStopCtx{2,1}(1,:),...
%     1:100, msRsdToPullStopStr{1,1}(1,:), msRsdToPullStopStr{2,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
% 
% figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(2,:), msRsdToPullStopCtx{2,1}(2,:),...
%     1:100, msRsdToPullStopStr{1,1}(2,:), msRsdToPullStopStr{2,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
% 
% figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(3,:), msRsdToPullStopCtx{2,1}(3,:),...
%     1:100, msRsdToPullStopStr{1,1}(3,:), msRsdToPullStopStr{2,1}(3,:), 'cmap', colorMap, 'transparency', 0.2);

%% Helper function
% Matrix interpolation function
    function [intMat] = intm(origMat, numbDataPoints)
        % 1-d interpolation of a matrix as specified by the number of data points
        % numbDataPoints = 100; % 20ms*100 = 2000ms
        x=1:size(origMat,2);
        xq=linspace(1,size(origMat,2),numbDataPoints);
        intMat = interp1(x,origMat',xq)';
        if size(intMat,1)>size(intMat,2)
            intMat = intMat'; % to ensure variable-by-time orientation
        end
    end
end
















