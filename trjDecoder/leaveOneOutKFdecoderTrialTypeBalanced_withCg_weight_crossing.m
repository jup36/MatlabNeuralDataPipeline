function  rez = leaveOneOutKFdecoderTrialTypeBalanced_withCg_weight_crossing(stateD_R, stateD_P, s)
%stateD = tmpState;

nd_R = {};  % container cell for all neural data
nd_P = {};

counter = 0;
if isfield(s.dat, 'spkCtxR')
    counter = counter + 1;
    nd_R{1, counter} = 'ctx';
    nd_R{2, counter} = s.dat.spkCtxR;
    nd_R{3, counter} = find(s.ctxI);
    valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd_R{2, counter}, 'un', 0)); % valid trials
    nd_R{4, counter} = mode(cellfun(@(a) size(a,1), nd_R{2, counter}(valTrI)));
    
    nd_P{1, counter} = 'ctx';
    nd_P{2, counter} = s.dat.spkCtxP;
    nd_P{3, counter} = find(s.ctxI);
    valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd_P{2, counter}, 'un', 0)); % valid trials
    nd_P{4, counter} = mode(cellfun(@(a) size(a,1), nd_P{2, counter}(valTrI)));
end

if isfield(s.dat, 'spkStrR')
    counter = counter + 1;
    nd_R{1, counter} = 'str';
    nd_R{2, counter} = s.dat.spkStrR;
    nd_R{3, counter} = find(s.strI);
    nd_R{4, counter} = mode(cellfun(@(a) size(a,1), nd_R{2, counter}(valTrI)));
    
    nd_P{1, counter} = 'str';
    nd_P{2, counter} = s.dat.spkStrP;
    nd_P{3, counter} = find(s.strI);
    nd_P{4, counter} = mode(cellfun(@(a) size(a,1), nd_P{2, counter}(valTrI)));
end

if isfield(s.dat, 'spkCgR')
    counter = counter + 1;
    nd_R{1, counter} = 'cg';
    nd_R{2, counter} = s.dat.spkCgR;
    if ~exist('valTrI', 'var')
        valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd_R{2, counter}, 'un', 0)); % valid trials
    end
    nd_R{4, counter} = mode(cellfun(@(a) size(a,1), nd_R{2, counter}(valTrI)));
    nd_R{3, counter} = (1:nd_R{4, counter})';
    
    nd_P{1, counter} = 'cg';
    nd_P{2, counter} = s.dat.spkCgP;
    nd_P{4, counter} = mode(cellfun(@(a) size(a,1), nd_P{2, counter}(valTrI)));
    nd_P{3, counter} = (1:nd_P{4, counter})';
end

if isfield(s.dat, 'laserIdx')
    laserI = s.dat.laserIdxR;
    stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, laserI, 'un', 0)); % stim trials
end

trainTrN = min(sum(valTrI & ~stmTrI))-1;

ctx_record = contains(cell2mat(nd_R(1,:)), 'ctx');
str_record = contains(cell2mat(nd_R(1,:)), 'str');
cg_record = contains(cell2mat(nd_R(1,:)), 'cg');

%% run
for r = 1:size(valTrI,1) % row: trials
    for c = 1:size(valTrI,2) % column: position/torque pairs
        trainI = true(size(valTrI,1),size(valTrI,2));
        if valTrI(r,c)
            trainI(r,c) = false; % to leave one trial out as a test trial
            testI = ~trainI;     % index for the one test trial left out
            valTrainI = trainI & valTrI & ~stmTrI; % valid train trial index (exclude stim trials)
            
            % get current train trials by resampling (to include the same # of trials for each trial type)
            trainTrs = cell(1,size(trainI,2));
            trainState_R = []; % current training data states (hand position: row 1-3, velocity: row 4-6
            trainState_P = []; % current training data states (hand position: row 1-3, velocity: row 4-6
            
            if ctx_record && str_record && cg_record  % Ctx, Str, Cg
                train_ctx_R = [];  % sample size matched between ctx and str
                train_str_R = [];
                train_cg_R = [];
                
                train_ctx_P = [];  % sample size matched between ctx and str
                train_str_P = [];
                train_cg_P = [];
            elseif ctx_record && str_record && ~cg_record  % Ctx, Str
                train_ctx_R = [];
                train_str_R = [];
                
                train_ctx_P = [];  % sample size matched between ctx and str
                train_str_P = [];
            elseif ~ctx_record && ~str_record && cg_record  % Cg only recording
                train_cg_R = [];
                train_cg_P = [];
            end
            
            for cc = 1:size(valTrainI,2)  % sample trials from each trial type
                tempValTr = find(valTrainI(:,cc));
                tempValTrRand = tempValTr(randperm(length(tempValTr)));
                trainTrs{1,cc} =tempValTrRand(1:trainTrN); %tempValTrRand(1:trainTrN); % take the set # of randomized trials from each trial type (column)
                
                trainState_R = [trainState_R; stateD_R(trainTrs{1,cc},cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                trainState_P = [trainState_P; stateD_P(trainTrs{1,cc},cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                
                train_obs_tt_R = cellfun(@(a) a(trainTrs{1,cc},cc), nd_R(2, :), 'un', 0); % select relevant trials reach
                train_obs_tt_P = cellfun(@(a) a(trainTrs{1,cc},cc), nd_P(2, :), 'un', 0); % select relevant trials pull
                
                % collect cell number matched observations: Reach
                if ctx_record && str_record && cg_record
                    train_ctx_R = [train_ctx_R; train_obs_tt_R{1}];  % ctx (matched between ctx & str)
                    train_str_R = [train_str_R; train_obs_tt_R{2}];  % str (matched between ctx & str)
                    train_cg_R = [train_cg_R; train_obs_tt_R{3}];  % cg
                    
                    train_ctx_P = [train_ctx_P; train_obs_tt_P{1}];  % ctx (matched between ctx & str)
                    train_str_P = [train_str_P; train_obs_tt_P{2}];  % str (matched between ctx & str)
                    train_cg_P = [train_cg_P; train_obs_tt_P{3}];  % cg
                    
                elseif ctx_record && str_record && ~cg_record
                    train_ctx_R = [train_ctx_R; train_obs_tt_R{1}];  % ctx (matched between ctx & str)
                    train_str_R = [train_str_R; train_obs_tt_R{2}];  % str (matched between ctx & str)
                    
                    train_ctx_P = [train_ctx_P; train_obs_tt_P{1}];  % ctx (matched between ctx & str)
                    train_str_P = [train_str_P; train_obs_tt_P{2}];  % str (matched between ctx & str)
                elseif ~ctx_record && ~str_record && cg_record
                    train_cg_R = [train_cg_R; train_obs_tt_R{1}];  % cg
                    train_cg_P = [train_cg_P; train_obs_tt_P{1}];  % cg
                end
            end
            clearvars cc
            
            if ctx_record && str_record && cg_record
                % CTX decoder: Pull phase / all CTX cells with matching weights 
                [rez.estStateCtxM{r,c}, rez.estStateCtxCov{r,c}, rez.estStateCtxM_mismatch{r,c}, rez.estStateCtxCov_mismatch{r,c}, rez.angle_weight_vectors_ctx{r,c}] = ...
                    train_test_KF_weight_switch(trainState_P, trainState_R, train_ctx_P, train_ctx_R, nd_P{2,1}, stateD_P{testI}, testI);          
                % STR decoder (# matched CTX-STR)
                [rez.estStateStrM{r,c}, rez.estStateStrCov{r,c}, rez.estStateStrM_mismatch{r,c}, rez.estStateStrCov_mismatch{r,c}, rez.angle_weight_vectors_str{r,c}] = ...
                    train_test_KF_weight_switch(trainState_P, trainState_R, train_str_P, train_str_R, nd_P{2,2}, stateD_P{testI}, testI);         
                % Cg decoder
                [rez.estStateCgM{r,c}, rez.estStateCgCov{r,c}, rez.estStateCgM_mismatch{r,c}, rez.estStateCgCov_mismatch{r,c}, rez.angle_weight_vectors_cg{r,c}] = ...
                    train_test_KF_weight_switch(trainState_P, trainState_R, train_cg_P, train_cg_R, nd_P{2,3}, stateD_P{testI}, testI);                 
            elseif ctx_record && str_record && ~cg_record
                % CTX decoder: Pull phase / all CTX cells with matching weights 
                [rez.estStateCtxM{r,c}, rez.estStateCtxCov{r,c}, rez.estStateCtxM_mismatch{r,c}, rez.estStateCtxCov_mismatch{r,c}, rez.angle_weight_vectors_ctx{r,c}] = ...
                    train_test_KF_weight_switch(trainState_P, trainState_R, train_ctx_P, train_ctx_R, nd_P{2,1}, stateD_P{testI}, testI);          
                % STR decoder (# matched CTX-STR)
                [rez.estStateStrM{r,c}, rez.estStateStrCov{r,c}, rez.estStateStrM_mismatch{r,c}, rez.estStateStrCov_mismatch{r,c}, rez.angle_weight_vectors_str{r,c}] = ...
                    train_test_KF_weight_switch(trainState_P, trainState_R, train_str_P, train_str_R, nd_P{2,2}, stateD_P{testI}, testI);   
            elseif ~ctx_record && ~str_record && cg_record
                % Cg decoder
                [rez.estStateCgM{r,c}, rez.estStateCgCov{r,c}, rez.estStateCgM_mismatch{r,c}, rez.estStateCgCov_mismatch{r,c}, rez.angle_weight_vectors_cg{r,c}] = ...
                    train_test_KF_weight_switch(trainState_P, trainState_R, train_cg_P, train_cg_R, nd_P{2,1}, stateD_P{testI}, testI);  
            end
        end
    end
end

clearvars r c

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estStateMean, estStateCov, estStateMean_mismatch, estStateCov_mismatch, angle_weight_vectors] = train_test_KF_weight_switch(train_stateP, train_stateR, train_obsP, train_obsR, observe, test_state, testIdx)

test_obs = observe{testIdx};

neuron_I = (1:size(test_obs,1))';  % cell Id

[NTrain,NTrType] = size(train_stateP);
numb_cell = mode(cell2mat(cellfun(@(a) size(a, 1), train_obsP, 'un', 0)));
%% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interStateM = 0;    % inter-state matrix
intraStateM = 0;    % intra-state matrix
obsStateM = 0;   % observed Ctx state matrix
obsStateM_mismatch = 0;
obsIntraStateM = 0; % observed intra-state matrix
count = 0;
for t = 1:size(train_stateP,1) % # of trial
    % for parameter A, state model describes how state evolves over time
    stateZ1 = train_stateP{t}(:,2:end); % Zt
    stateZ2 = train_stateP{t}(:,1:end-1); % Zt-1
    tmp1 = stateZ1*stateZ2'; % Zt*Zt-1' (6-by-6 matrix)
    interStateM = interStateM+tmp1; % sum Zt*Zt-1'
    tmp2 = stateZ2*stateZ2'; % Zt-1*Zt-1' (6-by-6 matrix)
    intraStateM = intraStateM+tmp2; % sum Zt-1*Zt-1'
    % for parameter C, observation model describes how observation relates to the state
    stateZ = train_stateP{t}; % state (position, velocity)
    stateZ_mismatch = train_stateR{t}; % state (position, velocity)
    obs = train_obsP{t}; % spike data (pull phase)
    obs_mismatch = train_obsR{t}; % spike data (reach phase)
    tmp3 = obs*stateZ'; % Xt*Zt' (Ctx cell-by-state mat)
    tmp3_mismatch = obs_mismatch*stateZ_mismatch'; 
    obsStateM = obsStateM+tmp3; % sum Ctx observation-state matrix
    obsStateM_mismatch = obsStateM_mismatch+tmp3_mismatch; 
    tmp4 = stateZ*stateZ'; % Zt*Zt'
    obsIntraStateM = obsIntraStateM+tmp4; % sum Zt*Zt'
    % for parameter Pi and V
    count=count+1;
    poolZStart(:,count) = train_stateP{t}(:,1);
end
clearvars t tt
A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
C=(obsStateM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx spikes, which defines how observation relates to the state
C_mismatch=(obsStateM_mismatch/count)/(obsIntraStateM/count);
Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
clearvars poolZStart

% Fit Q and R (variance)
sumQM = 0;
sumRM = 0;
sumRM_mismatch = 0;
count = 0;
for t = 1:size(train_stateP,1) % # of trial
    % for parameter Q, covariance of the state model
    stateZ1 = train_stateP{t}(:,2:end); % Zt
    stateZ2 = train_stateP{t}(:,1:end-1); % Zt-1
    tmp5 = stateZ1-A*stateZ2; % Zt-A*Zt-1
    tmp5 = tmp5*tmp5'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
    count = count + size(stateZ1,2);
    sumQM = sumQM + tmp5; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
    % for parameter R of Ctx spikes, covariance of the observation model
    stateZ = train_stateP{t}; % Zt
    stateZ_mismatch = train_stateR{t}; 
    obs = train_obsP{t}; % spike data Pull
    obs_mismatch = train_obsR{t}; % spike data Reach
    tmp6 = obs-C*stateZ; % Xt-C*Zt
    tmp6 = tmp6*tmp6'; % (Xt-C*Zt)(Xt-C*Zt)
    tmp6_mismatch = obs_mismatch-C_mismatch*stateZ_mismatch; 
    tmp6_mismatch = tmp6_mismatch*tmp6_mismatch'; % (Xt-C*Zt)(Xt-C*Zt)
    sumRM = sumRM+tmp6; % Sum(Xt-C*Zt)(Xt-C*Zt)
    sumRM_mismatch = sumRM_mismatch+tmp6_mismatch; 
end
clearvars t tt
Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
R = sumRM/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
R_mismatch = sumRM_mismatch/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
%% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curTrLength = size(test_state, 2); % position, velocity in 20ms bins
% Initialization
mu=Pi;   % sample mean
mu_mismatch=Pi; 
sigma=V; % sample covariance
sigma_mismatch = V; 
clearvars estStateMean estStateCov_*
valCellI_match = ~(sum(C,2)==0) & ~isnan(sum(C,2));
valCellI_mismatch = ~(sum(C_mismatch,2)==0) & ~isnan(sum(C_mismatch,2));
valCellI = valCellI_match & valCellI_mismatch; 

C_Val = C(valCellI,:); % to drop the cells with no spike at all from the mapping C, if not it leads to singular matrix warning
R_Val = R(valCellI,valCellI);

C_Val_mismatch = C_mismatch(valCellI,:); % to drop the cells with no spike at all from the mapping C, if not it leads to singular matrix warning
R_Val_mismatch = R_mismatch(valCellI,valCellI); 

angle_weight_vectors = angleTwoVectors(C_Val, C_Val_mismatch); 

for b = 1:curTrLength % # of timebins
    % One-step prediction by ctx and str data separately
    mu = A*mu; % A is the coefficient matrix for the state model that maps the previous states to current states
    sigma = A*sigma*A'+Q; % Q is the covariance of the state model
    
    mu_mismatch = A*mu_mismatch; 
    sigma_mismatch = A*sigma_mismatch*A'+Q; % Q is the covariance of the state model
    
    % compute the Kalman gain (needs to separately computed for ctx and str)
    K = sigma*C_Val'/(C_Val*sigma*C_Val'+R_Val); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model
    K_mismatch = sigma_mismatch*C_Val_mismatch'/(C_Val_mismatch*sigma_mismatch*C_Val_mismatch'+R_Val_mismatch); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model    
    
    % update the state by observation
    curObsX = test_obs(neuron_I,b); % take the current ctx spikes bin-by-bin
    curObsX = curObsX(valCellI,1);
    
    mu = mu + K*(curObsX-C_Val*mu);
    mu_mismatch = mu_mismatch + K_mismatch*(curObsX-C_Val_mismatch*mu_mismatch);
    
    sigma = sigma - K*C_Val*sigma;
    sigma_mismatch = sigma_mismatch - K_mismatch*C_Val_mismatch*sigma_mismatch; 
    
    estStateMean(:,b)=mu;
    estStateMean_mismatch(:,b) = mu_mismatch; 
    
    estStateCov(:,:,b)=sigma;
    estStateCov_mismatch(:,:,b)=sigma_mismatch;

end
clearvars b

end

function [representState,representStateCut,representIntState] = getRepresentativeStateEst(refState,estState,valTrId)
%refState = stateD;
%estState = estStateCtxMean;
%valTrId = valTrI;
nKv = mode(cell2mat(cellfun(@(a) size(a,1),refState(valTrId),'un',0)));
for rr = 1:size(refState,1)
    for cl = 1:size(refState,2)
        if valTrId(rr,cl)
            if sum(cell2mat(cellfun(@(a) sum(sum(isnan(a)))>0, estState(rr, cl, :), 'un', 0)))==0
                
                distToRef = squeeze(cell2mat(cellfun(@(a) sum(sum(sqrt((a-refState{rr,cl}).^2))), estState(rr,cl,:),'un',0))); % average across ctx resampled trials
                representState{rr,cl} = estState{rr,cl,find(distToRef==min(distToRef),1,'first')};
                tmpState = representState{rr,cl};
                tmpCutState = nan(nKv,50);
                if size(tmpState,2)>=5
                    tmpCutState(:,1:min(size(tmpState,2),50)) = tmpState(:,1:min(size(tmpState,2),50));
                    representStateCut{rr,cl} = tmpCutState;
                    representIntState{rr,cl} = intm(tmpState,50); % original estimated trajectory with interpolation
                end
            else
                representState{rr,cl} = NaN;
                representStateCut{rr,cl} = NaN;
                representIntState{rr,cl} = NaN;
            end
        end
    end
    clearvars r c
end
end


function c_indexed = cell_indexing(c_to_index, index)
c_indexed = cell(size(c_to_index));
val_c_I = cell2mat(cellfun(@(a) ~isempty(a), c_to_index, 'un', 0));
val_c_indexed = cellfun(@(a) a(index, :), c_to_index(val_c_I), 'un', 0);
c_indexed(val_c_I) = deal(val_c_indexed);
end

