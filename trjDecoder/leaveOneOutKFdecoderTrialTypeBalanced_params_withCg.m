function  rez = leaveOneOutKFdecoderTrialTypeBalanced_params_withCg(stateD, s, resample)
%stateD = tmpState;

nd = {};  % container cell for all neural data
counter = 0;
if isfield(s.dat, 'spkCtxR')
    counter = counter + 1;
    nd{1, counter} = 'ctx';
    nd{2, counter} = s.dat.spkCtxR;
    nd{3, counter} = find(s.ctxI);
    valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd{2, counter}, 'un', 0)); % valid trials
    nd{4, counter} = mode(cellfun(@(a) size(a,1), nd{2, counter}(valTrI)));
end

if isfield(s.dat, 'spkStrR')
    counter = counter + 1;
    nd{1, counter} = 'str';
    nd{2, counter} = s.dat.spkStrR;
    nd{3, counter} = find(s.strI);
    nd{4, counter} = mode(cellfun(@(a) size(a,1), nd{2, counter}(valTrI)));
end

if isfield(s.dat, 'spkCgR')
    counter = counter + 1;
    nd{1, counter} = 'cg';
    nd{2, counter} = s.dat.spkCgR;
    if ~exist('valTrI', 'var')
        valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd{2, counter}, 'un', 0)); % valid trials
    end
    nd{4, counter} = mode(cellfun(@(a) size(a,1), nd{2, counter}(valTrI)));
    nd{3, counter} = (1:nd{4, counter})';
end

if isfield(s.dat, 'laserIdx')
    laserI = s.dat.laserIdxR;
    stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, laserI, 'un', 0)); % stim trials
end

trainTrN = min(sum(valTrI&~stmTrI))-1;

ctx_record = contains(cell2mat(nd(1,:)), 'ctx');
str_record = contains(cell2mat(nd(1,:)), 'str');
cg_record = contains(cell2mat(nd(1,:)), 'cg');

if ctx_record && str_record && cg_record  % Ctx, Str, Cg
    nd_ctx_str = cellfun(@(a, b) [a; b], nd{2,1}, nd{2,2}, 'un', 0);
    nd_ctx_str_cg = cellfun(@(a, b, c) [a; b; c], nd{2,1}, nd{2,2}, nd{2,3}, 'un', 0);
elseif ctx_record && str_record && ~cg_record
    nd_ctx_str = cellfun(@(a, b) [a; b], nd{2,1}, nd{2,2}, 'un', 0);
end

%% run
for i = 1:resample  % repeat resampling trials
    Ns = cell2mat(nd(4, :)); % cell numbers
    
    rand_C = cellfun(@randperm, nd(4,:), 'un', 0);
    
    if ctx_record && str_record && cg_record  % Ctx, Str, Cg
        rs_n = min(Ns(1:2));  % number of cells to match between Ctx and Str
        rs_n_all = min(Ns);  % number of cells to match among Ctx, Str, and Cg
        rand_rs_n = cellfun(@(a) a(1:rs_n), rand_C(1:2), 'un', 0);
        rand_rs_n_all = cellfun(@(a) a(1:rs_n_all), rand_C, 'un', 0);
    elseif ctx_record && str_record && ~cg_record  % Ctx, Str
        rs_n = min(Ns(1:2));
        rand_rs_n = cellfun(@(a) a(1:rs_n), rand_C(1:2), 'un', 0);
    end
    
    %%
    for r = 1:size(valTrI,1) % row: trials
        for c = 1:size(valTrI,2) % column: position/torque pairs
            trainI = true(size(valTrI,1),size(valTrI,2));
            if valTrI(r,c)
                trainI(r,c) = false; % to leave one trial out as a test trial
                testI = ~trainI;     % index for the one test trial left out
                valTrainI = trainI & valTrI & ~stmTrI; % valid train trial index (exclude stim trials)
                
                % get current train trials by resampling (to include the same # of trials for each trial type)
                trainTrs = cell(1,size(trainI,2));
                trainState = []; % current training data states (hand position: row 1-3, velocity: row 4-6
                
                if ctx_record && str_record && cg_record  % Ctx, Str, Cg
                    train_ctx = [];  % sample size matched between ctx and str
                    train_str = [];
                    train_cg = [];
                    train_ctx1 = [];  % sample size matched to cg
                    train_str1 = [];
                    train_ctx_str = [];  % ctx, str combined
                    train_ctx_str_cg = [];  % ctx, str, cg combined
                    
                elseif ctx_record && str_record && ~cg_record  % Ctx, Str
                    train_ctx = [];
                    train_str = [];
                    train_ctx_str = [];
                    
                elseif ~ctx_record && ~str_record && cg_record  % Cg only recording
                    train_cg = [];
                end
                
                trainNd = cell(1, size(nd,2));
                
                for cc = 1:size(valTrainI,2)  % sample trials from each trial type
                    tempValTr = find(valTrainI(:,cc));
                    tempValTrRand = tempValTr(randperm(length(tempValTr)));
                    trainTrs{1,cc} =tempValTrRand(1:trainTrN); %tempValTrRand(1:trainTrN); % take the set # of randomized trials from each trial type (column)
                    
                    trainState = [trainState; stateD(trainTrs{1,cc},cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                    
                    train_obs_tt = cellfun(@(a) a(trainTrs{1,cc},cc), nd(2, :), 'un', 0); % select relevant trials
                    
                    % collect cell number matched observations
                    if ctx_record && str_record && cg_record
                        train_ctx = [train_ctx; cellfun(@(a) a(rand_rs_n{1}, :), train_obs_tt{1}, 'un', 0)];  % ctx (matched between ctx & str)
                        train_str = [train_str; cellfun(@(a) a(rand_rs_n{2}, :), train_obs_tt{2}, 'un', 0)];  % str (matched between ctx & str)
                        train_ctx1 = [train_ctx1; cellfun(@(a) a(rand_rs_n_all{1}, :), train_obs_tt{1}, 'un', 0)];  % ctx (matched to cg)
                        train_str1 = [train_str1; cellfun(@(a) a(rand_rs_n_all{2}, :), train_obs_tt{2}, 'un', 0)];  % str (matched to cg)
                        train_cg = [train_cg; train_obs_tt{3}];  % cg
                        train_ctx_str = [train_ctx_str; cellfun(@(a, b) [a; b], train_obs_tt{1}, train_obs_tt{2}, 'un', 0)];  % ctx, str combined
                        train_ctx_str_cg = [train_ctx_str_cg; cellfun(@(a, b, c) [a; b; c], train_obs_tt{1}, train_obs_tt{2}, train_obs_tt{3}, 'un', 0)];  % ctx, str, cg combined
                    elseif ctx_record && str_record && ~cg_record
                        train_ctx = [train_ctx; cellfun(@(a) a(rand_rs_n{1}, :), train_obs_tt{1}, 'un', 0)];  % ctx (matched between ctx & str)
                        train_str = [train_str; cellfun(@(a) a(rand_rs_n{2}, :), train_obs_tt{2}, 'un', 0)];  % str (matched between ctx & str)
                        train_ctx_str = [train_ctx_str; cellfun(@(a, b) [a; b], train_obs_tt{1}, train_obs_tt{2}, 'un', 0)];  % ctx, str combined
                        %train_ctx_str = [train_ctx_str; cellfun(@(a, b) [a; b], train_ctx, train_str, 'un', 0)];  % ctx, str combined
                    elseif ~ctx_record && ~str_record && cg_record
                        train_cg = [train_cg; train_obs_tt{3}];  % cg
                    end
                end
                clearvars cc
                
                if ctx_record && str_record && cg_record
                    % CTX decoder (# matched CTX-STR)
                    [rez.estStateCtxM{r,c,i}, rez.estStateCtxCov{r,c,i}, rez.coeff_Ctx{r,c,i}] = train_test_KF(trainState, train_ctx, nd{2,1}, stateD{testI}, testI, ...
                        rand_rs_n{1}, nd{3,1}(rand_rs_n{1}), length(s.ctxI)); % ctx
                    % STR decoder (# matched CTX-STR)
                    [rez.estStateStrM{r,c,i}, rez.estStateStrCov{r,c,i}, rez.coeff_Str{r,c,i}] = train_test_KF(trainState, train_str, nd{2,2}, stateD{testI}, testI, ...
                        rand_rs_n{2}, nd{3,2}(rand_rs_n{2}), length(s.strI)); % str
                    % CTX decoder (# matched with Cg)
                    [rez.estStateCtxM_cg{r,c,i}, rez.estStateCtxCov_cg{r,c,i}, ~] = train_test_KF(trainState, train_ctx1, nd{2,1}, stateD{testI}, testI, ...
                        rand_rs_n_all{1}, nd{3,1}(rand_rs_n_all{1}), length(s.ctxI)); % ctx
                    % STR decoder (# matched with Cg)
                    [rez.estStateStrM_cg{r,c,i}, rez.estStateStrCov_cg{r,c,i}, ~] = train_test_KF(trainState, train_str1, nd{2,2}, stateD{testI}, testI, ...
                        rand_rs_n_all{2}, nd{3,2}(rand_rs_n_all{2}), length(s.strI)); % str
                    % Cg decoder
                    [rez.estStateCgM{r,c,i}, rez.estStateCgCov{r,c,i}, ~] = train_test_KF(trainState, train_cg, nd{2,3}, stateD{testI}, testI, ...
                        nd{3,3}, nd{3,3}, length(nd{3,3})); % Cg
                    % CTX-STR decoder
                    [rez.estStateCtxStrM{r,c,i}, rez.estStateCtxStrCov{r,c,i}, rez.coeff_CtxStr{r,c,i}] = train_test_KF(trainState, train_ctx_str, nd_ctx_str, stateD{testI}, testI, ...
                        cell2mat(nd(3,1:2)'), cell2mat(nd(3,1:2)'), length(s.ctxI)); % Ctx & Str
                    % CTX-STR-Cg decoder
                    [rez.estStateCtxStrCgM{r,c,i}, rez.estStateCtxStrCgCov{r,c,i}, ~] = train_test_KF(trainState, train_ctx_str_cg, nd_ctx_str_cg, stateD{testI}, testI); % Ctx & Str & Cg
                elseif ctx_record && str_record && ~cg_record
                    % CTX decoder (# matched CTX-STR)
                    [rez.estStateCtxM{r,c,i}, rez.estStateCtxCov{r,c,i}, rez.coeff_Ctx{r,c,i}] = train_test_KF(trainState, train_ctx, nd{2,1}, stateD{testI}, testI, ...
                        rand_rs_n{1}, nd{3,1}(rand_rs_n{1}), length(s.ctxI)); % ctx
                    % STR decoder (# matched CTX-STR)
                    [rez.estStateStrM{r,c,i}, rez.estStateStrCov{r,c,i}, rez.coeff_Str{r,c,i}] = train_test_KF(trainState, train_str, nd{2,2}, stateD{testI}, testI, ...
                        rand_rs_n{2}, nd{3,2}(rand_rs_n{2}), length(s.strI)); % str
                    % CTX-STR decoder
                    [rez.estStateCtxStrM{r,c,i}, rez.estStateCtxStrCov{r,c,i}, rez.coeff_CtxStr{r,c,i}] = train_test_KF(trainState, train_ctx_str, nd_ctx_str, stateD{testI}, testI, ...
                        cell2mat(nd(3,1:2)'), cell2mat(nd(3,1:2)'), length(s.ctxI)); % Ctx & Str
                elseif ~ctx_record && ~str_record && cg_record
                    % Cg decoder
                    [rez.estStateCgM{r,c,i}, rez.estStateCgCov{r,c,i}, ~] = train_test_KF(trainState, train_cg, nd{2,3}, stateD{testI}, testI, ...
                        nd{3,3}, nd{3,3}, length(nd{3,3})); % Cg
                end
            end
        end
    end
    fprintf('finished resampling iteration #%d\n', i)
end
clearvars r c


if ctx_record && str_record && cg_record  % Ctx, Str, Cg
    [rez.rpr_est_ctx, rez.rpr_est_cut_ctx, rez.rpr_est_intp_ctx] = getRepresentativeStateEst(stateD, rez.estStateCtxM, valTrI);
    [rez.rpr_est_str, rez.rpr_est_cut_str, rez.rpr_est_intp_str] = getRepresentativeStateEst(stateD, rez.estStateStrM, valTrI);
    [rez.rpr_est_ctx_str, rez.rpr_est_cut_ctx_str, rez.rpr_est_intp_ctx_str] = getRepresentativeStateEst(stateD, rez.estStateCtxStrM, valTrI);
    [rez.rpr_est_ctx_str_cg, rez.rpr_est_cut_ctx_str_cg, rez.rpr_est_intp_ctx_str_cg] = getRepresentativeStateEst(stateD, rez.estStateCtxStrCgM, valTrI);
    [rez.rpr_est_cg, rez.rpr_est_cut_cg, rez.rpr_est_intp_cg] = getRepresentativeStateEst(stateD, rez.estStateCgM, valTrI);
elseif ctx_record && str_record && ~cg_record  % Ctx, Str
    [rez.rpr_est_ctx, rez.rpr_est_cut_ctx, rez.rpr_est_intp_ctx] = getRepresentativeStateEst(stateD, rez.estStateCtxM, valTrI);
    [rez.rpr_est_str, rez.rpr_est_cut_str, rez.rpr_est_intp_str] = getRepresentativeStateEst(stateD, rez.estStateStrM, valTrI);
    [rez.rpr_est_ctx_str, rez.rpr_est_cut_ctx_str, rez.rpr_est_intp_ctx_str] = getRepresentativeStateEst(stateD, rez.estStateCtxStrM, valTrI);
elseif ~ctx_record && ~str_record && cg_record  % Cg only recording
    [rez.rpr_est_cg, rez.rpr_est_cut_cg, rez.rpr_est_intp_cg] = getRepresentativeStateEst(stateD, rez.estStateCgM, valTrI);
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estStateMean, estStateCov, C_Val_out] = train_test_KF(train_state, train_obs, observe, test_state, testIdx, varargin)
% train_state: cell array containing trial-by-trial state variables
% train_obs: cell array containing trial-by-trial observation variables (neural data)
% obs: cell array containing all observation variable (neural data)
% test_state: array containing held-out state variable
% testIdx: index indicating the held out trial
% cI: resampled cell id (e.g., rand_rs_n{1});
% ref_cI: resampled cells' coordinates in the original binSpkCountCell;

test_obs = observe{testIdx};

if nargin == 5
   neuron_I = (1:size(test_obs,1))';  % cell Id
elseif nargin > 5
   cI = varargin{1};
   neuron_I = cI';  % used cell Id
   neuron_I_ref = varargin{2};  % the location of used cells in the reference vector
   size_ref_vec = varargin{3}; 
end

[NTrain,NTrType] = size(train_state);
numb_cell = mode(cell2mat(cellfun(@(a) size(a, 1), train_obs, 'un', 0)));
%% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interStateM = 0;    % inter-state matrix
intraStateM = 0;    % intra-state matrix
obsStateM = 0;   % observed Ctx state matrix
obsIntraStateM = 0; % observed intra-state matrix
count = 0;
for t = 1:size(train_state,1) % # of trial
    % for parameter A, state model describes how state evolves over time
    stateZ1 = train_state{t}(:,2:end); % Zt
    stateZ2 = train_state{t}(:,1:end-1); % Zt-1
    tmp1 = stateZ1*stateZ2'; % Zt*Zt-1' (6-by-6 matrix)
    interStateM = interStateM+tmp1; % sum Zt*Zt-1'
    tmp2 = stateZ2*stateZ2'; % Zt-1*Zt-1' (6-by-6 matrix)
    intraStateM = intraStateM+tmp2; % sum Zt-1*Zt-1'
    % for parameter C, observation model describes how observation relates to the state
    stateZ = train_state{t}; % state (position, velocity)
    obs = train_obs{t}; % spike data
    tmp3 = obs*stateZ'; % Xt*Zt' (Ctx cell-by-state mat)
    obsStateM = obsStateM+tmp3; % sum Ctx observation-state matrix
    tmp4 = stateZ*stateZ'; % Zt*Zt'
    obsIntraStateM = obsIntraStateM+tmp4; % sum Zt*Zt'
    % for parameter Pi and V
    count=count+1;
    poolZStart(:,count) = train_state{t}(:,1);
end
clearvars t tt
A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
C=(obsStateM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx spikes, which defines how observation relates to the state
Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
clearvars poolZStart

% Fit Q and R (variance)
sumQM = 0;
sumRM = 0;
count = 0;
for t = 1:size(train_state,1) % # of trial
    % for parameter Q, covariance of the state model
    stateZ1 = train_state{t}(:,2:end); % Zt
    stateZ2 = train_state{t}(:,1:end-1); % Zt-1
    tmp5 = stateZ1-A*stateZ2; % Zt-A*Zt-1
    tmp5 = tmp5*tmp5'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
    count = count + size(stateZ1,2);
    sumQM = sumQM + tmp5; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
    % for parameter R of Ctx spikes, covariance of the observation model
    stateZ = train_state{t}; % Zt
    obs = train_obs{t}; % spike data
    tmp6 = obs-C*stateZ; % Xt-C*Zt
    tmp6 = tmp6*tmp6'; % (Xt-C*Zt)(Xt-C*Zt)
    sumRM = sumRM+tmp6; % Sum(Xt-C*Zt)(Xt-C*Zt)
end
clearvars t tt
Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
R = sumRM/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
%% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curTrLength = length(test_state); % position, velocity in 20ms bins
% Initialization
mu=Pi;   % sample mean
sigma=V; % sample covariance
clearvars estStateMean estStateCov_*
valCellI = ~(sum(C,2)==0) & ~isnan(sum(C,2));

C_Val = C(valCellI,:); % to drop the cells with no spike at all from the mapping C, if not it leads to singular matrix warning
R_Val = R(valCellI,valCellI);

if nargin == 5 % no need to further reference
    tempC = nan(numb_cell, size(C_Val,2));
    tempC(valCellI, :) = C_Val;
    C_Val_out = tempC;
elseif nargin > 5 % use reference to place coefficients correspondingly
    tempC = nan(size_ref_vec, size(C_Val,2));
    tempC(neuron_I_ref(valCellI),:) = C_Val;
    C_Val_out = tempC;
end

for b = 1:curTrLength % # of timebins
    % One-step prediction by ctx and str data separately
    mu = A*mu; % A is the coefficient matrix for the state model that maps the previous states to current states
    sigma = A*sigma*A'+Q; % Q is the covariance of the state model
    
    % compute the Kalman gain (needs to separately computed for ctx and str)
    K = sigma*C_Val'/(C_Val*sigma*C_Val'+R_Val); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model
    
    % update the state by observation
    curObsX = test_obs(neuron_I,b); % take the current ctx spikes bin-by-bin
    curObsX = curObsX(valCellI,1);
    
    mu = mu + K*(curObsX-C_Val*mu);
    sigma = sigma - K*C_Val*sigma;
    estStateMean(:,b)=mu;
    estStateCov(:,:,b)=sigma;
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
            distToRef = squeeze(cell2mat(cellfun(@(a) sum(sum(sqrt((a-refState{rr,cl}).^2))),estState(rr,cl,:),'un',0))); % average across ctx resampled trials
            representState{rr,cl} = estState{rr,cl,find(distToRef==min(distToRef),1,'first')};
            tmpState = representState{rr,cl};
            tmpCutState = nan(nKv,50);
            if size(tmpState,2)>=5
                tmpCutState(:,1:min(size(tmpState,2),50)) = tmpState(:,1:min(size(tmpState,2),50));
                representStateCut{rr,cl} = tmpCutState;
                representIntState{rr,cl} = intm(tmpState,50); % original estimated trajectory with interpolation
            end
        end
    end
end
clearvars r c
end

function c_indexed = cell_indexing(c_to_index, index)
c_indexed = cell(size(c_to_index));
val_c_I = cell2mat(cellfun(@(a) ~isempty(a), c_to_index, 'un', 0));
val_c_indexed = cellfun(@(a) a(index, :), c_to_index(val_c_I), 'un', 0);
c_indexed(val_c_I) = deal(val_c_indexed);
end

