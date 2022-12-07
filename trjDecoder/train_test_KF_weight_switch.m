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