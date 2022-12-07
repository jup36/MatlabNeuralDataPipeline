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
curTrLength = size(test_state, 2); % position, velocity in 20ms bins
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