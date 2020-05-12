fileDirectory = '/Users/parkj/Documents/MATLAB/MatlabNeuralDataPipeline/Neural_signal_processing_BY/PS9_Kalman_filter';
cd(fileDirectory)
load('ps9_data.mat')

[NTrain, NDirect]=size(train_trial); % # of trials (91), # of reach directions (8)
NTest=size(test_trial,1); % # of test trials

%% %%%%%%%%%%%%%%%%%%%% Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% get the 1) neuron-by-timebin spike count matrix, 2) binned state matrix containing position & velocity 
timeStep=20; % 20 ms bins

% Train dataset
for trialIX=1:NTrain % # of trials along each direction: 91
    for directIX=1:NDirect % # of directions: 8
        curSpikes=train_trial(trialIX,directIX).spikes; % current spikes
        curPos=train_trial(trialIX,directIX).handPos;   % current hand positions
        timeWinE=timeStep:timeStep:size(curSpikes,2);   % timebins 
        curSpikes=curSpikes(:,1:timeWinE(end-1))';      % timebin x neuron spike matrix
        curPos=curPos(:,1:timeWinE(end));               % current hand positions truncated
        % Take spike counts in each bin
        spikeMtr=reshape(curSpikes,timeStep,[],size(curSpikes,2));    % get spike count matrices (time x timesteps x neurons)
        train(trialIX,directIX).spikeCount=squeeze(sum(spikeMtr,1))'; % get the neuron-by-timebin (get the binned spike counts by sum) spike count matrix
        % Compute a 4-dimensional arm state
        pos=curPos(1:2,timeWinE);                       % get the binned hand position
        vel=(pos(:,2:end)-pos(:,1:end-1))/timeStep;     % get the binned velocity
        train(trialIX,directIX).state=[pos(:,1:end-1); vel]; % get the state matrix (position & velocity)
    end
end

% Test dataset
for trialIX=1:NTest % # of trials along each direction: 91
    for directIX=1:NDirect % # of directions: 8
        curSpikes=test_trial(trialIX,directIX).spikes; % current spikes
        curPos=test_trial(trialIX,directIX).handPos;   % current hand positions
        timeWinE=timeStep:timeStep:size(curSpikes,2);  % timebins
        curSpikes=curSpikes(:,1:timeWinE(end-1))';     % timebin x neuron spike matrix
        curPos=curPos(:,1:timeWinE(end));              % current hand positions truncated
        % Take spike counts in each bin
        spikeMtr=reshape(curSpikes,timeStep,[],size(curSpikes,2));      % get spike count matrices (time x timesteps x neurons)
        test(trialIX,directIX).spikeCount=squeeze(sum(spikeMtr,1))';    % get the neuron-by-timebin spike count matrix
        % Compute a 4-dimensional arm state
        pos=curPos(1:2,timeWinE);                           % get the binned hand position
        vel=(pos(:,2:end)-pos(:,1:end-1))/timeStep;         % get the binned velocity
        test(trialIX,directIX).state=[pos(:,1:end-1); vel]; % get the 4-d state matrix (position & velocity)
    end
end
clear train_trial;
clear test_trial;

%% %%%%%%%%%%%%%%%%%%%% Training Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Fit A, Pi, V and C
interStateMtr=0;    % inter-state matrix
intraStateMtr=0;    % intra-state matrix
obsStateMtr=0;      % observed state matrix
obsIntraStateMtr=0; % observed intrastate matrix
count=0;
countPool=0;
for trialIX=1:NTrain % # of Trial
    for directIX=1:NDirect % # of Direction
        % for parameter A, state model describes how state evolves over time
        stateZ1=train(trialIX,directIX).state(:,2:end);     % Zt
        stateZ2=train(trialIX,directIX).state(:,1:end-1);   % Zt-1
        tmp1=stateZ1*stateZ2';                              % Zt*Zt-1' (4-by-4 matrix)
        interStateMtr=interStateMtr+tmp1;                   % Sum Zt*Zt-1'
        tmp2=stateZ2*stateZ2';                              % Zt-1*Zt-1' (4-by-4 matrix)
        intraStateMtr=intraStateMtr+tmp2;                   % Sum Zt-1*Zt-1'
        % for parameter C, observation model describes how observation relates to the state
        stateZ=train(trialIX,directIX).state;       % State (position, velocity) 
        obserX=train(trialIX,directIX).spikeCount;  % Observation (spike data)
        tmp3=obserX*stateZ';                    % Xt*Zt' (97-by-4 matrix)
        obsStateMtr=obsStateMtr+tmp3;           % Sum Xt*Zt'
        tmp4=stateZ*stateZ';                    % Zt*Zt'
        obsIntraStateMtr=obsIntraStateMtr+tmp4; % Sum Zt*Zt'
        % for parameter Pi and V
        count=count+1;
        poolZStart(:,count)=train(trialIX,directIX).state(:,1);
    end
end
A=(interStateMtr/count)*inv(intraStateMtr/count);   % Slope for the state model, which defines how state evolves over time
C=(obsStateMtr/count)*inv(obsIntraStateMtr/count);  % Slope for the observation model, which defines how observation relates to the state
Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
V=cov(poolZStart');    % Covariance for the initial state (sample covariance)

% Fit Q and R
sumQMtr=0;
sumRMtr=0; 
count=0;
for trialIX=1:NTrain % # of Trial
    for directIX=1:NDirect % # of Direction
        % for parameter Q, covariance of the state model 
        stateZ1=train(trialIX,directIX).state(:,2:end);
        stateZ2=train(trialIX,directIX).state(:,1:end-1);
        tmp5=stateZ1-A*stateZ2; % Zt-A*Zt-1
        tmp5=tmp5*tmp5';        % (Zt-A*Zt-1)(Zt-A*Zt-1)'
        count=count+size(stateZ1,2);
        sumQMtr=sumQMtr+tmp5;   % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
        % for parameter R, covariance of the observation model
        stateZ=train(trialIX,directIX).state;       % Zt
        obserX=train(trialIX,directIX).spikeCount;  % Xt (Neuron-by-timebins)
        tmp6=obserX-C*stateZ;   % Xt-C*Zt 
        tmp6=tmp6*tmp6';        % (Xt-C*Zt)(Xt-C*Zt)'
        sumRMtr=sumRMtr+tmp6;    % Sum (Xt-C*Zt)(Xt-C*Zt)'   
    end
end
Q=sumQMtr/(count);                % count 1/(t-1), count should correspond to 2 to T
R=sumRMtr/(count+NTrain*NDirect); % count 1/t, count should correspond to 1 to T

%% %%%%%%%%%%%%%%%%%%%%%%%% Test Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for trialIX=1:NTest % # of trials along each direction: 91
    for directIX=1:NDirect % # of directions: 8
        curTrialLen=size(test(trialIX,directIX).state,2); % position, velocity in 20 ms bins
        % Initialization
        mu=Pi;   % sample mean
        sigma=V; % sample covariance
        clear curEstStateMean; 
        clear curEstStateCov;
        for timeIX=1:curTrialLen % # of current timebins (varies trial-by-trial)
            % One-step Prediction
            mu=A*mu;            % A is the coefficient matrix for the state model that maps the previous states to current states
            sigma=A*sigma*A'+Q; % Q is the covariance of the state model
            % compute the Kalman gain
            K=sigma*C'*inv(C*sigma*C'+R); % Kalman gain
            % update the state
            curObservX=test(trialIX,directIX).spikeCount(:,timeIX); % take the current spike-count bin
            mu=mu+K*(curObservX-C*mu);
            sigma=sigma-K*C*sigma;
            
            curEstStateMean(:,timeIX)=mu;   % current estimation (bin-by-bin) of the state (mean)
            curEstStateCov(:,:,timeIX)=sigma; % current estimation (bin-by-bin) of the state (covariance)
        end
        test(trialIX,directIX).estStateMean=curEstStateMean;
        test(trialIX,directIX).estStateCov=curEstStateCov;
    end
end

%% %%%%%%%%%%%%%%%%%%%% Plot the selected trials %%%%%%%%%%%%%%%%%%%%%%%%%%
selTrialList={{1,1};{1,2};{1,3};{1,4};{1,5};{1,6};{1,7};{1,8}};
for trialIX=1:size(selTrialList,1) % # of selected trials
    selTrial=selTrialList{trialIX};
    curTrial=test(selTrial{1},selTrial{2}); % (Trial,Direction)
    figure(trialIX);
    titleStr=sprintf('Trial#: %d   Direction#: %d',selTrial{1},selTrial{2});
    
    % 2-d trajectory x and y
    subplot(1,3,1);
    plot(curTrial.state(1,:),curTrial.state(2,:),'-ok','LineWidth',2); % actual state trajectory
    hold on;
    plot(curTrial.estStateMean(1,:),curTrial.estStateMean(2,:),'-or','LineWidth',2); % estimated state trajectory
    % draw covariance ellipses 
    for timeIX=1:length(curTrial.estStateMean(1,:)) % # of bins
        hold on;
        func_plotEllipse(curTrial.estStateMean(1:2,timeIX),curTrial.estStateCov(1:2,1:2,timeIX));
    end
    xlabel('Horz-Pos'); ylabel('Vert-Pos'); legend('Org','Est');
    title(titleStr);
    
    % x (horizontal) trajectory
    subplot(1,3,2); % horizontal position
    selState=1;     % horizontal position
    curStd=(squeeze(sqrt(curTrial.estStateCov(selState,selState,:))))'; % standard deviation (Horz-Pos)
    curTimeIX=20*(1:length(curStd)); % time steps
    
    plot(curTimeIX,curTrial.state(selState,:),'-*k');                        % actual x-traj
    hold on; plot(curTimeIX,curTrial.estStateMean(selState,:),'-*r');        % estimated x-traj
    hold on; plot(curTimeIX,curTrial.estStateMean(selState,:)+curStd,'--r'); % estimated x-traj + std
    hold on; plot(curTimeIX,curTrial.estStateMean(selState,:)-curStd,'--r'); % estimated x-traj - std
    xlabel('time'); ylabel('Horz-Pos'); legend('Org','Est');
    title(titleStr);
    
    % y (vertical) trajectory
    selState=2; % vertical position
    subplot(1,3,3); % vertical position
    selState=2;     % vertical position
    curStd=(squeeze(sqrt(curTrial.estStateCov(selState,selState,:))))'; % standard deviation (Horz-Pos)
    curTimeIX=20*(1:length(curStd)); % time steps
    
    plot(curTimeIX,curTrial.state(selState,:),'-*k');                        % actual y-traj
    hold on; plot(curTimeIX,curTrial.estStateMean(selState,:),'-*r');        % estimated y-traj
    hold on; plot(curTimeIX,curTrial.estStateMean(selState,:)+curStd,'--r'); % estimated y-traj + std
    hold on; plot(curTimeIX,curTrial.estStateMean(selState,:)-curStd,'--r'); % estimated y-traj - std
    xlabel('time'); ylabel('Vert-Pos'); legend('Org','Est');
    title(titleStr);
end












