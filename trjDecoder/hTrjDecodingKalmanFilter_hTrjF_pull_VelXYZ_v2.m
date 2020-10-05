function hTrjDecodingKalmanFilter_hTrjF_pull_VelXYZ_v2(filePath, saveName, plotlogic)
%This decodes kinematics of mouse 3-d hand movement trajectories (X,Y,Z)
% using cross-validated (leave-a-trial-out) Kalman filter decoding.
% Note these changes in the version 2: 
% 1) preprocessing code: preprocessForKFdecoder_HTrjF_rotate_reachPull.m 
% 2) Using trajectories on the rotated coordinates (X, Y are now purely X,Y)
% 3) Also decoding from all neurons in both strital and cortical populations
% 4) Using all valid trials available, so the # of trials per trial type
%    might differ.  

%plot logic = 0;
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)

kfDir = dir('preprocessKFdecodeHTrjCtxStr_reachPull_rotate_baseSub*');
%kfDir = dir('preprocessKFdecodeHTrjCtxStr_reachPull_rotate*');
load(fullfile(kfDir.folder,kfDir.name),'s')
%fileName = 'preprocessKFdecodeHTrjCtxStr_WR40_081919.mat';
%load(fullfile(filePath,fileName),'s')

figSaveDir = fullfile(filePath,'Figure','KalmanFilter_decoding_pull');
if ~(isfolder(fullfile(filePath,'Figure','KalmanFilter_decoding_pull')))
    mkdir(fullfile(filePath,'Figure','KalmanFilter_decoding_pull'))
end

resample = 20; %100;
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxP, 'un', 0));
stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, s.dat.laserIdx, 'un', 0)); % stim trials
%trainTrN = min(sum(valTrI & ~stmTrI))-1;

ctxCnumb = max(unique(cell2mat(cellfun(@(a) size(a,1), s.dat.spkCtxP, 'un', 0)))); % # of cortex cells
strCnumb = max(unique(cell2mat(cellfun(@(a) size(a,1), s.dat.spkStrP, 'un', 0)))); % # of striatal cells
minCnumb = min(ctxCnumb, strCnumb); % # of cells to be included

%% select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
for r = 1:size(s.dat.stateP,1)
    for c = 1:size(s.dat.stateP,2)
        if ~isempty(s.dat.stateP{r,c})
            s.dat.stateP{r,c} = s.dat.stateP{r,c}(4:6,:); % X, Y, Z velocity (row 1:3 for position, row 4:6 for velocity)
        end
    end
end
clearvars r c

%%%%%%%%%%%%%%% NOTE THAT s.dat.stateP contains hand trajectories already
%%%%%%%%%%%%%%% median subtracted and rotated 

%% leave-a-trial-out decoding using Kalman Filter (heavy-lifting part)
for i = 1:resample %resample % repeat resampling trials
    randCtxI = randperm(ctxCnumb);
    randStrI = randperm(strCnumb);
    ctxI = randCtxI(1:minCnumb); % ctx cells for this iteration
    strI = randStrI(1:minCnumb); % str cells for this iteration
    
    for r = 1:size(s.dat.spkCtxP,1) % row: trials
        for c = 1:size(s.dat.spkCtxP,2) % column: position/torque pairs
            trainI = true(size(valTrI,1),size(valTrI,2));
            if valTrI(r,c)
                trainI(r,c) = false; % to leave one trial out as a test trial
                testI = ~trainI;     % index for the one test trial left out
                valTrainI = trainI & valTrI & ~stmTrI; % valid train trial index (exclude stim trials)
                
                % get current train trials by resampling (to include the same # of trials for each trial type)
                currTrainTrs = cell(1,size(trainI,2)); % zeros(trainTrN,size(trainI,2));
                currTrainState = []; % current training data states (hand position: row 1-3, velocity: row 4-6)
                currTrainCtx = [];
                currTrainStr = [];
                currTrainCtxStr = []; % decoding using both populations
                for cc = 1:size(valTrainI,2) % sample trials from each trial type
                    currTrainTrs{1,cc} = find(valTrainI(:,cc)); % use all available trials
                    currTrainState = [currTrainState; s.dat.stateP(currTrainTrs{1,cc},cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                    currTrainCtx = [currTrainCtx; cellfun(@(a) a(ctxI,:), s.dat.spkCtxP(currTrainTrs{1,cc},cc), 'un', 0)]; % ctx spike mat with cells resampled to match # of cells
                    currTrainStr = [currTrainStr; cellfun(@(a) a(strI,:), s.dat.spkStrP(currTrainTrs{1,cc},cc), 'un', 0)]; % str spike mat with cells resampled to match # of cells
                    currTrainCtxStr = [currTrainCtxStr; cellfun(@(a,b) [a(ctxI,:); b(strI,:)], s.dat.spkCtxP(currTrainTrs{1,cc},cc), s.dat.spkStrP(currTrainTrs{1,cc},cc), 'un', 0)]; 
                end
                clearvars cc
                [NTrain,NTrType] = size(currTrainState);
                
                %% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                interStateM = 0;    % inter-state matrix
                intraStateM = 0;    % intra-state matrix
                obsStateCtxM = 0;   % observed Ctx state matrix
                obsStateStrM = 0;   % observed Str state matrix
                obsStateCtxStrM = 0;% observed CtxStr state matrix 
                obsIntraStateM = 0; % observed intra-state matrix
                count = 0;
                countPool = 0;
                for t = 1:size(currTrainState,1) % # of trial
                    % for parameter A, state model describes how state evolves over time
                    stateZ1 = currTrainState{t}(:,2:end); % Zt
                    stateZ2 = currTrainState{t}(:,1:end-1); % Zt-1
                    tmp1 = stateZ1*stateZ2'; % Zt*Zt-1' (6-by-6 matrix)
                    interStateM = interStateM+tmp1; % sum Zt*Zt-1'
                    tmp2 = stateZ2*stateZ2'; % Zt-1*Zt-1' (6-by-6 matrix)
                    intraStateM = intraStateM+tmp2; % sum Zt-1*Zt-1'
                    % for parameter C, observation model describes how observation relates to the state
                    stateZ = currTrainState{t}; % state (position, velocity)
                    obsCtx = currTrainCtx{t}; % ctx spike data
                    obsStr = currTrainStr{t}; % str spike data
                    obsCtxStr = currTrainCtxStr{t}; % ctx and str spike data 
                    tmp3 = obsCtx*stateZ'; % Xt*Zt' (Ctx cell-by-state mat)
                    obsStateCtxM = obsStateCtxM+tmp3; % sum Ctx observation-state matrix
                    tmp4 = obsStr*stateZ'; % Xt*Zt' (Str cell-by-state mat)
                    obsStateStrM = obsStateStrM+tmp4; % sum Str observation-state matrix
                    tmp9 = obsCtxStr*stateZ'; % Xt*Zt' (CtxStr cell-by-state mat)
                    obsStateCtxStrM = obsStateCtxStrM+tmp9; % sum CtxStr observation-state matrix
                    tmp5 = stateZ*stateZ'; % Zt*Zt'
                    obsIntraStateM = obsIntraStateM+tmp5; % sum Zt*Zt'
                    % for parameter Pi and V
                    count=count+1;
                    poolZStart(:,count) = currTrainState{t}(:,1);
                end
                clearvars t tt
                A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
                C_ctx=(obsStateCtxM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx spikes, which defines how observation relates to the state
                C_str=(obsStateStrM/count)/(obsIntraStateM/count); % Slope for the observation model for str spikes
                C_ctxstr = (obsStateCtxStrM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx str spikes 
                Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
                V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
                clearvars poolZStart
                
                % Fit Q and R (variance)
                sumQM = 0;
                sumRM_ctx = 0;
                sumRM_str = 0;
                sumRM_ctxstr = 0; 
                count = 0;
                for t = 1:size(currTrainState,1) % # of trial
                    % for parameter Q, covariance of the state model
                    stateZ1 = currTrainState{t}(:,2:end); % Zt
                    stateZ2 = currTrainState{t}(:,1:end-1); % Zt-1
                    tmp6 = stateZ1-A*stateZ2; % Zt-A*Zt-1
                    tmp6 = tmp6*tmp6'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
                    count = count + size(stateZ1,2);
                    sumQM = sumQM + tmp6; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
                    % for parameter R of Ctx spikes, covariance of the observation model
                    stateZ = currTrainState{t}; % Zt
                    obsCtx = currTrainCtx{t}; % ctx spike data
                    tmp7 = obsCtx-C_ctx*stateZ; % Xt_ctx-C_ctx*Zt
                    tmp7 = tmp7*tmp7'; % (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                    sumRM_ctx = sumRM_ctx+tmp7; % Sum (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                    % for parameter R of Str spikes, covariance of the observation model
                    obsStr = currTrainStr{t}; % str spike data
                    tmp8 = obsStr-C_str*stateZ; % Xt_str-C_str*Zt
                    tmp8 = tmp8*tmp8'; % (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    sumRM_str = sumRM_str+tmp8; % Sum (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    % for parameter R of CtxStr spikes, covariance of the observation model
                    obsCtxStr = currTrainCtxStr{t}; % ctx str spike data
                    tmp10 = obsCtxStr-C_ctxstr*stateZ; % Xt_ctxstr-C_ctxstr*Zt
                    tmp10 = tmp10*tmp10'; % (Xt_ctxstr-C_ctxstr*Zt)(Xt_ctxstr-C_ctxstr*Zt)
                    sumRM_ctxstr = sumRM_ctxstr+tmp10; % Sum (Xt_ctxstr-C_ctxstr*Zt)(Xt_ctxstr-C_ctxstr*Zt)    
                end
                clearvars t tt
                Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
                R_ctx = sumRM_ctx/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                R_str = sumRM_str/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                R_ctxstr = sumRM_ctxstr/(count+NTrain*NTrType); 
                %% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                curTrLength = size(s.dat.stateP{testI},2); % position, velocity in 20ms bins
                % Initialization
                mu_ctx=Pi;   % sample mean
                sigma_ctx=V; % sample covariance
                mu_str=Pi;   % sample mean
                sigma_str=V; % sample covariance
                mu_ctxstr = Pi; % sample mean
                sigma_ctxstr = V; % sample covariance
                
                clearvars curEstStateMean*;
                clearvars curEstStateCov*;
                valCellICtx = ~(sum(C_ctx,2)==0)&~isnan(sum(C_ctx,2));
                valCellIStr = ~(sum(C_str,2)==0)&~isnan(sum(C_str,2));
                valCellICtxStr = ~(sum(C_ctxstr,2)==0)&~isnan(sum(C_ctxstr,2));
                
                C_ctxVal = C_ctx(valCellICtx,:); % to drop the cells with no spike at all from the mapping C, if not it leads to singular matrix warning
                R_ctxVal = R_ctx(valCellICtx,valCellICtx); 
                C_strVal = C_str(valCellIStr,:); 
                R_strVal = R_str(valCellIStr,valCellIStr); 
                C_ctxstrVal = C_ctxstr(valCellICtxStr,:);
                R_ctxstrVal = R_ctxstr(valCellICtxStr,valCellICtxStr);
                
                for b = 1:curTrLength % # of timebins
                    % One-step prediction by ctx and str data separately
                    mu_ctx = A*mu_ctx; % A is the coefficient matrix for the state model that maps the previous states to current states
                    sigma_ctx = A*sigma_ctx*A'+Q; % Q is the covariance of the state model
                    mu_str = A*mu_str;
                    sigma_str = A*sigma_str*A'+Q;
                    mu_ctxstr = A*mu_ctxstr; 
                    sigma_ctxstr = A*sigma_ctxstr*A'+Q;
                    
                    % compute the Kalman gain (needs to separately computed for ctx and str)
                    K_ctx = sigma_ctx*C_ctxVal'/(C_ctxVal*sigma_ctx*C_ctxVal'+R_ctxVal); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model
                    K_str = sigma_str*C_strVal'/(C_strVal*sigma_str*C_strVal'+R_strVal); % *inv(C_str*sigma*C_str'+R_str); % Kalman gain for str observation model
                    K_ctxstr = sigma_ctxstr*C_ctxstrVal'/(C_ctxstrVal*sigma_ctxstr*C_ctxstrVal'+R_ctxstrVal); 
                    
                    % update the state by ctx observation
                    curObsX_ctx = s.dat.spkCtxP{testI}(ctxI,b); % take the current ctx spikes bin-by-bin
                    curObsX_ctx = curObsX_ctx(valCellICtx,1); 
                    
                    mu_ctx = mu_ctx + K_ctx*(curObsX_ctx-C_ctxVal*mu_ctx);
                    sigma_ctx = sigma_ctx - K_ctx*C_ctxVal*sigma_ctx;
                    curEstStateMean_ctx(:,b)=mu_ctx;
                    curEstStateCov_ctx(:,:,b)=sigma_ctx;
                    
                    % update the state by str observation
                    curObsX_str = s.dat.spkStrP{testI}(strI,b); % take the current str spikes bin-by-bin
                    curObsX_str = curObsX_str(valCellIStr,1); 
                    
                    mu_str = mu_str + K_str*(curObsX_str-C_strVal*mu_str);
                    sigma_str = sigma_str - K_str*C_strVal*sigma_str;
                    curEstStateMean_str(:,b)=mu_str;
                    curEstStateCov_str(:,:,b)=sigma_str;
                    
                    % update the state by ctxstr observation
                    curObsX_ctxstr = [s.dat.spkCtxP{testI}(ctxI,b);s.dat.spkStrP{testI}(strI,b)]; % take the current ctx spikes bin-by-bin
                    curObsX_ctxstr = curObsX_ctxstr([valCellICtx;valCellIStr],1); 
                    
                    mu_ctxstr = mu_ctxstr + K_ctxstr*(curObsX_ctxstr-C_ctxstrVal*mu_ctxstr);
                    sigma_ctxstr = sigma_ctxstr - K_ctxstr*C_ctxstrVal*sigma_ctxstr;
                    curEstStateMean_ctxstr(:,b)=mu_ctxstr;
                    curEstStateCov_ctxstr(:,:,b)=sigma_ctxstr;
                    
                end
                s.dat.estStateCtxMean{r,c,i} = curEstStateMean_ctx;
                s.dat.estStateCtxCov{r,c,i} = curEstStateCov_ctx;
                
                s.dat.estStateStrMean{r,c,i} = curEstStateMean_str;
                s.dat.estStateStrCov{r,c,i} = curEstStateCov_str;
                
                s.dat.estStateCtxStrMean{r,c,i} = curEstStateMean_ctxstr;
                s.dat.estStateCtxStrCov{r,c,i} = curEstStateCov_ctxstr;
                
            end
        end
    end
    fprintf('finished iteration# %d\n', i)
end
clearvars r c

%% get the # of kinematic variables
valCellI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.stateP, 'un', 0));
sizeCell = cellfun(@size, s.dat.stateP, 'un', 0);
valSizeC = cell2mat(sizeCell(valCellI));
nKv = unique(valSizeC(:,1)); % the # of kinematic variables

%% take average across estimated trajectories of resampling, interpolate to match trajectory lengths
for rr = 1:size(s.dat.estStateCtxMean,1) % trials (row)
    for cc = 1:size(s.dat.estStateCtxMean,2) % trial-types (column)
        % cortex estimated state
        s.dat.statePCtx{rr,cc} = nanmedian(cell2mat(s.dat.estStateCtxMean(rr,cc,:)),3); % average across ctx resampled trials
        tmpCtx = s.dat.statePCtx{rr,cc};
        stateCtx{rr,cc} = nan(nKv,100); % original estimated trajectory without interpolation
        if size(tmpCtx,2)>=5
            stateCtx{rr,cc}(:,1:min(size(tmpCtx,2),100)) = tmpCtx(:,1:min(size(tmpCtx,2),100));
            stateCtxInt{rr,cc} = intm(tmpCtx,100); % original estimated trajectory with interpolation
        end
        
        % striatum estimated state
        s.dat.statePStr{rr,cc} = nanmedian(cell2mat(s.dat.estStateStrMean(rr,cc,:)),3); % average across str resampled trials
        tmpStr = s.dat.statePStr{rr,cc};
        stateStr{rr,cc} = nan(nKv,100); % original estimated trajectory without interpolation
        if size(tmpStr,2)>=5
            stateStr{rr,cc}(:,1:min(size(tmpStr,2),100)) = tmpStr(:,1:min(size(tmpStr,2),100));
            stateStrInt{rr,cc} = intm(tmpStr,100); % original estimated trajectory with interpolation
        end
        
        % cortexstriatum estimated state
        s.dat.statePCtxStr{rr,cc} = nanmedian(cell2mat(s.dat.estStateCtxStrMean(rr,cc,:)),3); % average across ctxstr resampled trials
        tmpCtxStr = s.dat.statePCtxStr{rr,cc};
        stateCtxStr{rr,cc} = nan(nKv,100); % original estimated trajectory without interpolation
        if size(tmpCtxStr,2)>=5
            stateCtxStr{rr,cc}(:,1:min(size(tmpCtxStr,2),100)) = tmpCtxStr(:,1:min(size(tmpCtxStr,2),100));
            stateCtxStrInt{rr,cc} = intm(tmpCtxStr,100); % original estimated trajectory with interpolation
        end
        
        % actual state
        tmpAct = s.dat.stateP{rr,cc};
        state{rr,cc} = nan(nKv,100);
        
        if size(tmpAct,2)>=5
            state{rr,cc}(:,1:min(size(tmpAct,2),100)) = tmpAct(:,1:min(size(tmpAct,2),100));
            stateInt{rr,cc} = intm(tmpAct,100);
        end
        
        % distance between the actual trajectory and the mean estimated cortex and striatum trajectory
        if ~isempty(s.dat.stateP{rr,cc})
            s.dat.dTrjCtx{rr,cc} = sqrt((s.dat.stateP{rr,cc}-s.dat.statePCtx{rr,cc}).^2);
            s.dat.dTrjStr{rr,cc} = sqrt((s.dat.stateP{rr,cc}-s.dat.statePStr{rr,cc}).^2);
            s.dat.dTrjStrStr{rr,cc} = sqrt((s.dat.stateP{rr,cc}-s.dat.statePCtxStr{rr,cc}).^2);
        else
            s.dat.dTrjCtx{rr,cc} = [];
            s.dat.dTrjStr{rr,cc} = [];
            s.dat.dTrjCtxStr{rr,cc} = [];
        end
    end
end
clearvars rr cc

%% Get the mean and sem trajectories for each direction-torque combinations
valCellI = cell2mat(cellfun(@(a) ~isempty(a), stateCtxInt, 'un', 0));
sizeCell = cellfun(@size, stateCtxInt, 'un', 0);
valSizeC = cell2mat(sizeCell(valCellI));
nTb = unique(valSizeC(:,2)); % the # of time bins

for cc = 1:size(s.dat.estStateCtxMean,2) % trial-types
    %% cortex-estimated interpolated trajectories
    % cortex estimated trajectory with interpolation NO STIM/LASER trials whole
    tmpKvTimeTrialCtxInt = cell2mat(reshape(stateCtxInt(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxInt{1,cc} = nanmean(tmpKvTimeTrialCtxInt,3);
    s.dat.trSCtxInt{1,cc} = nanstd(tmpKvTimeTrialCtxInt,0,3)./sqrt(size(tmpKvTimeTrialCtxInt,3));
    
    % cortex estimated trajectory with interpolation STIM/LASER trials
    tmpKvTimeTrialCtxIntLaser = cell2mat(reshape(stateCtxInt(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxIntLaser{1,cc} = nanmean(tmpKvTimeTrialCtxIntLaser,3);
    s.dat.trSCtxIntLaser{1,cc} = nanstd(tmpKvTimeTrialCtxIntLaser,0,3)./sqrt(size(tmpKvTimeTrialCtxIntLaser,3));
    
    %% striatum-estimated interpolated trajectories
    % striatum estimated trajectory with interpolation NO STIM/LASER trials whole
    tmpKvTimeTrialStrInt = cell2mat(reshape(stateStrInt(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMStrInt{1,cc} = nanmean(tmpKvTimeTrialStrInt,3);
    s.dat.trSStrInt{1,cc} = nanstd(tmpKvTimeTrialStrInt,0,3)./sqrt(size(tmpKvTimeTrialStrInt,3));
    
    % striatum estimated trajectory with interpolation STIM/LASER trials
    tmpKvTimeTrialStrIntLaser = cell2mat(reshape(stateStrInt(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMStrIntLaser{1,cc} = nanmean(tmpKvTimeTrialStrIntLaser,3);
    s.dat.trSStrIntLaser{1,cc} = nanstd(tmpKvTimeTrialStrIntLaser,0,3)./sqrt(size(tmpKvTimeTrialStrIntLaser,3));
    
    %% cortexstriatum-estimated interpolated trajectories
    % cortexstriatum estimated trajectory with interpolation NO STIM/LASER trials whole
    tmpKvTimeTrialCtxStrInt = cell2mat(reshape(stateCtxStrInt(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxStrInt{1,cc} = nanmean(tmpKvTimeTrialCtxStrInt,3);
    s.dat.trSCtxStrInt{1,cc} = nanstd(tmpKvTimeTrialCtxStrInt,0,3)./sqrt(size(tmpKvTimeTrialCtxStrInt,3));
    
    % cortexstriatum estimated trajectory with interpolation STIM/LASER trials
    tmpKvTimeTrialCtxStrIntLaser = cell2mat(reshape(stateCtxStrInt(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxStrIntLaser{1,cc} = nanmean(tmpKvTimeTrialCtxStrIntLaser,3);
    s.dat.trSCtxStrIntLaser{1,cc} = nanstd(tmpKvTimeTrialCtxStrIntLaser,0,3)./sqrt(size(tmpKvTimeTrialCtxStrIntLaser,3));
    
    %% actual interpolated trajectories
    % actual hand trajectory with interpolation NO STIM/LASER trials
    tmpKvTimeTrialActInt = cell2mat(reshape(stateInt(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMActInt{1,cc} = nanmean(tmpKvTimeTrialActInt,3);
    s.dat.trSActInt{1,cc} = nanstd(tmpKvTimeTrialActInt,0,3)./sqrt(size(tmpKvTimeTrialActInt,3));
    
    % actual hand trajectory with interpolation STIM/LASER trials
    tmpKvTimeTrialActIntLaser = cell2mat(reshape(stateInt(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMActIntLaser{1,cc} = nanmean(tmpKvTimeTrialActIntLaser,3);
    s.dat.trSActIntLaser{1,cc} = nanstd(tmpKvTimeTrialActIntLaser,0,3)./sqrt(size(tmpKvTimeTrialActIntLaser,3));
    
    %% uninterpolated trajectories cortex striatum actual
    %% cortex estimated trajectory without interpolation NO STIM/LASER trials
    tmpKvTimeTrialCtx = cell2mat(reshape(stateCtx(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtx{1,cc} = nanmean(tmpKvTimeTrialCtx,3);
    s.dat.trSCtx{1,cc} = nanstd(tmpKvTimeTrialCtx,0,3)./sqrt(size(tmpKvTimeTrialCtx,3));
    % cortex estimated trajectory without interpolation STIM/LASER trials
    tmpKvTimeTrialCtxLaser = cell2mat(reshape(stateCtx(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxLaser{1,cc} = nanmean(tmpKvTimeTrialCtxLaser,3);
    s.dat.trSCtxLaser{1,cc} = nanstd(tmpKvTimeTrialCtxLaser,0,3)./sqrt(size(tmpKvTimeTrialCtxLaser,3));
    
    %% striatum estimated trajectory without interpolation NO STIM/LASER trials
    tmpKvTimeTrialStr = cell2mat(reshape(stateStr(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMStr{1,cc} = nanmean(tmpKvTimeTrialStr,3);
    s.dat.trSStr{1,cc} = nanstd(tmpKvTimeTrialStr,0,3)./sqrt(size(tmpKvTimeTrialStr,3));
    % striatum estimated trajectory without interpolation STIM/LASER trials
    tmpKvTimeTrialStrLaser = cell2mat(reshape(stateStr(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMStrLaser{1,cc} = nanmean(tmpKvTimeTrialStrLaser,3);
    s.dat.trSStrLaser{1,cc} = nanstd(tmpKvTimeTrialStrLaser,0,3)./sqrt(size(tmpKvTimeTrialStrLaser,3));
    
    %% cortexstriatum estimated trajectory without interpolation NO STIM/LASER trials
    tmpKvTimeTrialCtxStr = cell2mat(reshape(stateCtxStr(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxStr{1,cc} = nanmean(tmpKvTimeTrialCtxStr,3);
    s.dat.trSCtxStr{1,cc} = nanstd(tmpKvTimeTrialCtxStr,0,3)./sqrt(size(tmpKvTimeTrialCtxStr,3));
    % cortexstriatum estimated trajectory without interpolation STIM/LASER trials
    tmpKvTimeTrialCtxStrLaser = cell2mat(reshape(stateCtxStr(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMCtxStrLaser{1,cc} = nanmean(tmpKvTimeTrialCtxStrLaser,3);
    s.dat.trSCtxStrLaser{1,cc} = nanstd(tmpKvTimeTrialCtxStrLaser,0,3)./sqrt(size(tmpKvTimeTrialCtxStrLaser,3));
    
    % actual hand trajectory without interpolation NO STIM/LASER trials
    tmpKvTimeTrialAct = cell2mat(reshape(state(valCellI(:,cc)&~stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMAct{1,cc} = nanmean(tmpKvTimeTrialAct,3);
    s.dat.trSAct{1,cc} = nanstd(tmpKvTimeTrialAct,0,3)./sqrt(size(tmpKvTimeTrialAct,3));
    % actual hand trajectory without interpolation STIM/LASER trials
    tmpKvTimeTrialActLaser = cell2mat(reshape(state(valCellI(:,cc)&stmTrI(:,cc),cc),1,1,[]));
    s.dat.trMActLaser{1,cc} = nanmean(tmpKvTimeTrialActLaser,3);
    s.dat.trSActLaser{1,cc} = nanstd(tmpKvTimeTrialActLaser,0,3)./sqrt(size(tmpKvTimeTrialActLaser,3));
end
clearvars cc

%% evaluate decodting with correlation and r-squred 
[corrRez.ctx,r2Rez.ctx] = trjDecodeEvalByTrialType(s.dat.stateP, s.dat.statePCtx); 
[corrRez.str,r2Rez.str] = trjDecodeEvalByTrialType(s.dat.stateP, s.dat.statePStr); 
[corrRez.ctxstr,r2Rez.ctxstr] = trjDecodeEvalByTrialType(s.dat.stateP, s.dat.statePCtxStr); 

%% save results
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrVel_pull_v2',saveName)),'s','corrRez','r2Rez') % last saved after training without stim trials 5/27 Wed 9pm
%load(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrVel_',saveName)),'s') % last saved after training without stim trials 5/27 Wed 9pm
%trjMovie([stateCtxCC_sm(:,2), stateStrCC_sm(:,2), stateCC_sm(:,2)]', figSaveDir, 'kfDecode_Ypos_CtxStrAct')

%% plot X, Y, Z trajectories separately for high versus low torques
if plotlogic == 1
colorMap = [[100 149 237]./255; [50 205 50]./255; [50 50 50]./255]; % colorMap for cortex and striatum
%% X trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,1}(1,:),s.dat.trSCtxInt{1,1}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,1}(1,:),s.dat.trSStrInt{1,1}(1,:), ...
        1:nTb,s.dat.trMActInt{1,1}(1,:),s.dat.trSActInt{1,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_whole_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,1}(1,:),s.dat.trSCtxIntLaser{1,1}(1,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,1}(1,:),s.dat.trSStrIntLaser{1,1}(1,:), ...
        1:nTb,s.dat.trMActIntLaser{1,1}(1,:),s.dat.trSActIntLaser{1,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_Laser_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% X trj, low torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,3}(1,:),s.dat.trSCtxInt{1,3}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,3}(1,:),s.dat.trSStrInt{1,3}(1,:), ...
        1:nTb,s.dat.trMActInt{1,3}(1,:),s.dat.trSActInt{1,3}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_whole_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    % X trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,3}(1,:),s.dat.trSCtxIntLaser{1,3}(1,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,3}(1,:),s.dat.trSStrIntLaser{1,3}(1,:), ...
        1:nTb,s.dat.trMActIntLaser{1,3}(1,:),s.dat.trSActIntLaser{1,3}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_Laser_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% X trj, high torque, position 1(left), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,2}(1,:),s.dat.trSCtxInt{1,2}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,2}(1,:),s.dat.trSStrInt{1,2}(1,:), ...
        1:nTb,s.dat.trMActInt{1,2}(1,:),s.dat.trSActInt{1,2}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_whole_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % X trj, high torque, position 1(left), Cortex, Striatum, Actual LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,2}(1,:),s.dat.trSCtxIntLaser{1,2}(1,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,2}(1,:),s.dat.trSStrIntLaser{1,2}(1,:), ...
        1:nTb,s.dat.trMActIntLaser{1,2}(1,:),s.dat.trSActIntLaser{1,2}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_Laser_HtP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    %% X trj, high torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,4}(1,:),s.dat.trSCtxInt{1,4}(1,:), ...
        1:nTb,s.dat.trMStrInt{1,4}(1,:),s.dat.trSStrInt{1,4}(1,:), ...
        1:nTb,s.dat.trMActInt{1,4}(1,:),s.dat.trSActInt{1,4}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_whole_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % X trj, high torque, position 2(right), Cortex, Striatum, Actual LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,4}(1,:),s.dat.trSCtxIntLaser{1,4}(1,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,4}(1,:),s.dat.trSStrIntLaser{1,4}(1,:), ...
        1:nTb,s.dat.trMActIntLaser{1,4}(1,:),s.dat.trSActIntLaser{1,4}(1,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('xVel_interp_Laser_HtP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Y trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,1}(2,:),s.dat.trSCtxInt{1,1}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,1}(2,:),s.dat.trSStrInt{1,1}(2,:), ...
        1:nTb,s.dat.trMActInt{1,1}(2,:),s.dat.trSActInt{1,1}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_whole_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,1}(2,:),s.dat.trSCtxIntLaser{1,1}(2,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,1}(2,:),s.dat.trSStrIntLaser{1,1}(2,:), ...
        1:nTb,s.dat.trMActIntLaser{1,1}(2,:),s.dat.trSActIntLaser{1,1}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_Laser_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Y trj, low torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,3}(2,:),s.dat.trSCtxInt{1,3}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,3}(2,:),s.dat.trSStrInt{1,3}(2,:), ...
        1:nTb,s.dat.trMActInt{1,3}(2,:),s.dat.trSActInt{1,3}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_whole_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, low torque, position 2(right), Cortex, Striatum, Actual whole trajectory LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,3}(2,:),s.dat.trSCtxIntLaser{1,3}(2,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,3}(2,:),s.dat.trSStrIntLaser{1,3}(2,:), ...
        1:nTb,s.dat.trMActIntLaser{1,3}(2,:),s.dat.trSActIntLaser{1,3}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_Laser_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Y trj, high torque, position 1(left), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,2}(2,:),s.dat.trSCtxInt{1,2}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,2}(2,:),s.dat.trSStrInt{1,2}(2,:), ...
        1:nTb,s.dat.trMActInt{1,2}(2,:),s.dat.trSActInt{1,2}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_whole_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, high torque, position 1(left), Cortex, Striatum, Actual whole trajectory LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,2}(2,:),s.dat.trSCtxIntLaser{1,2}(2,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,2}(2,:),s.dat.trSStrIntLaser{1,2}(2,:), ...
        1:nTb,s.dat.trMActIntLaser{1,2}(2,:),s.dat.trSActIntLaser{1,2}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_Laser_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Y trj, high torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,4}(2,:),s.dat.trSCtxInt{1,4}(2,:), ...
        1:nTb,s.dat.trMStrInt{1,4}(2,:),s.dat.trSStrInt{1,4}(2,:), ...
        1:nTb,s.dat.trMActInt{1,4}(2,:),s.dat.trSActInt{1,4}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_whole_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
    
    % Y trj, high torque, position 2(right), Cortex, Striatum, Actual whole trajectory LASER
    figure;
    boundedline(1:nTb,s.dat.trMCtxIntLaser{1,4}(2,:),s.dat.trSCtxIntLaser{1,4}(2,:), ...
        1:nTb,s.dat.trMStrIntLaser{1,4}(2,:),s.dat.trSStrIntLaser{1,4}(2,:), ...
        1:nTb,s.dat.trMActIntLaser{1,4}(2,:),s.dat.trSActIntLaser{1,4}(2,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-12 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('yVel_interp_Laser_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Z trj, low torque, position 1(left), Cortex, Striatum, Actual whole trajectory
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,1}(3,:),s.dat.trSCtxInt{1,1}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,1}(3,:),s.dat.trSStrInt{1,1}(3,:), ...
        1:nTb,s.dat.trMActInt{1,1}(3,:),s.dat.trSActInt{1,1}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('zVel_interp_whole_ltP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Z trj, low torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,3}(3,:),s.dat.trSCtxInt{1,3}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,3}(3,:),s.dat.trSStrInt{1,3}(3,:), ...
        1:nTb,s.dat.trMActInt{1,3}(3,:),s.dat.trSActInt{1,3}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('zVel_interp_whole_ltP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Z trj, high torque, position 1(left), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,2}(3,:),s.dat.trSCtxInt{1,2}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,2}(3,:),s.dat.trSStrInt{1,2}(3,:), ...
        1:nTb,s.dat.trMActInt{1,2}(3,:),s.dat.trSActInt{1,2}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('zVel_interp_whole_htP1',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% Z trj, high torque, position 2(right), Cortex, Striatum, Actual
    figure;
    boundedline(1:nTb,s.dat.trMCtxInt{1,4}(3,:),s.dat.trSCtxInt{1,4}(3,:), ...
        1:nTb,s.dat.trMStrInt{1,4}(3,:),s.dat.trSStrInt{1,4}(3,:), ...
        1:nTb,s.dat.trMActInt{1,4}(3,:),s.dat.trSActInt{1,4}(3,:), 'cmap', colorMap, 'transparency', 0.2);
    ylim([-6 8])
    print(fullfile(filePath,'Figure','KalmanFilter_decoding_pull',strcat('zVel_interp_whole_htP2',saveName)),'-dpdf','-painters','-bestfit')
    close;
        
    %% plot concatenated
    timeX = 0;
    figure; hold on;
    for c = 1:size(s.dat.estStateCtxMean,2)
        for r = 1:20 % just to include the first block only per trial type
            if ~isempty(s.dat.stateP{r,c})
                ctxTrjX = s.dat.statePCtx{r,c}(1,:); % Ctx X trj (left-right, horizontal hand position)
                strTrjX = s.dat.statePStr{r,c}(1,:); % Str X trj
                actTrjX = s.dat.stateP{r,c}(1,:); % actual trj
                
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    function [corrRez,r2Rez] = trjDecodeEvalByTrialType(state, stateEst)
        %Computes pearson correlation and r-squared metrics between actual and
        % estimated states, each input should be trial-by-trialType cell, whose cell
        % element should be variables-by-timeBin.
        % state = s.dat.stateP; % actual state
        % stateEst = s.dat.statePCtx; % estimated state
        s1 = cell2mat(reshape(state,[],1)')'; % state mat
        s2 = cell2mat(reshape(stateEst,[],1)')'; % stateEst mat
        
        r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
        r2Rez.all = r2(s1,s2); % r^2 across all trials
        
        corrRez.all = corr(s1,s2,'Rows','complete'); % corr without smoothing
        corrRez.all_sm = corr(smooth2a(s1,3,0),smooth2a(s2,3,0),'Rows','complete'); % corr with smoothing across time bins (smoothing seems to hurt usually)
        
        for tt = 1:size(state,2) % trial types
            s1tt = cell2mat(reshape(state(:,tt),[],1)')';
            s2tt = cell2mat(reshape(stateEst(:,tt),[],1)')';
            tempCorr = corr(s1tt,s2tt,'Rows','complete');
            tempCorr_sm = corr(smooth2a(s1tt,3,0),smooth2a(s2tt,3,0),'Rows','complete');
            tempR2 = r2(s1tt,s2tt);
            
            switch tt
                case 1
                    corrRez.lelt = tempCorr;
                    corrRez.lelt_sm = tempCorr_sm;
                    r2Rez.lelt = tempR2;
                case 2
                    corrRez.leht = tempCorr;
                    corrRez.leht_sm = tempCorr_sm;
                    r2Rez.leht = tempR2;
                case 3
                    corrRez.rilt = tempCorr;
                    corrRez.rilt_sm = tempCorr_sm;
                    r2Rez.rilt = tempR2;
                case 4
                    corrRez.riht = tempCorr;
                    corrRez.riht_sm = tempCorr_sm;
                    r2Rez.riht = tempR2;
            end
        end
    end

end
















