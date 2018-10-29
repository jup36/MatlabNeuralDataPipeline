%This script implements Kalman filter that maps neuronal spike activity to 
% movement kinematic states - specifically joystick position (not x, y trajectories 
% but more like summed absolute position of the joystick, i.e., amplitude) and velocity. 
% k-fold cross validation is used (crossvaltrials.m). 
% use dms data only

clear all; clear functions; clc;

fileDirectory = '/Volumes/RAID2/parkj/NeuralData/WR23_acc_dms_090117/Matfiles';
cd(fileDirectory)

addpath(genpath('/Volumes/RAID2/parkj/MATLAB'))
% addpath(genpath('/Users/parkj/Desktop/SpikeGLX-master 2/MATLAB-SDK')) % the folder with spikeGLX matlab scripts  
% addpath(genpath('/Volumes/RAID2/parkj/MATLAB/TONIC-master')) 

%binName =  'WR23_090117_Ldms_g0_t0.imec.ap.bin'; % .bin filename (e.g. 'JP_WR16_2130um_dura_062117_g0_t0.nidq.bin')
%meta    = ReadMeta(binName, '/Volumes/RAID2/parkj/NeuralData/WR23_Ldms_090117/RawData');         % get the meta data (structure)

%load('binSpkCountCTXreachStartWR23_090117', 'binSpkCountCTX1ms') % load the ctx binned spike counts
load('binSpkCountSTRreachStartWR23_090117', 'binSpkCountSTR1ms') % load the str binned spike counts
load('BehVariables','xpos1','ypos1', 'pos1', 'vel1', 'ts')       % load behavioral data - x & y positions aligned to the reachStart
binSC = binSpkCountSTR1ms;      % binned spike count structure

%% Variables for preprocessing
uv.neuralTime = -1000:3000-1;    % time points of the neural spike trains
uv.behavTime  = -200:1:1500-1;   % time points of the behavioral trajectories 
uv.kfTime     = [-100 1200];     % time points to be included in KF fit
uv.numbTrial  = size(ypos1,1);   % the number of trials based on the behavioral data 
uv.timeStep   = 20;              % timesteps or binsize in ms (stick with 20 ms as a default) 
uv.fold       = 5;               % for k-fold crossvalidation
uv.frLowCut   = 1;               % to exclude units with too sparse firing from the dataset  
uv.neuralTimeIdx = uv.neuralTime>=uv.kfTime(1) & uv.neuralTime<uv.kfTime(2); % time index for spike train data
uv.behavTimeIdx  = uv.behavTime>=uv.kfTime(1)  & uv.behavTime<uv.kfTime(2);  % time index for behavioral trajectory data

%% Construct the train and test datasets for 4-fold cross-validation 
% Rearrange the entire dataset
uCount = 0;
stackBinSC = zeros(uv.numbTrial,sum(uv.neuralTimeIdx),1); % stacked binned spike count mat

% This loop selects and stack units with all trials valid
for u = 1:length(binSC) % increment units
    if size(binSC(u).SpkCountMat,1)==uv.numbTrial % in case unit has got all trials
        tmpSCmtr = full(binSC(u).SpkCountMat(:,uv.neuralTimeIdx));
        tmpFR    = (sum(tmpSCmtr(:))*1000)/(size(tmpSCmtr,1)*size(tmpSCmtr,2)); % current unit's FR
        if tmpFR >= uv.frLowCut % in case FR greater or equal to 1 Hz
            uCount=uCount+1; % count
            stackBinSC(:,:,uCount)=full(binSC(u).SpkCountMat(:,uv.neuralTimeIdx)); % stack the unit's trial-by-trial psth
        end
    end
end
clearvars uCount binSpkCount*

% set trials for train vs. test data sets for k-fold cross-validation 
[ testTrialLogic ] = crossvaltrials( size(stackBinSC,1), uv.fold ); % get the train trial index for each cross-validation fold
testTrials         = cell2mat(testTrialLogic(:,2)); % testTrials
% clearvars -except binSC stackBinSC uv testTrialLogic xpos1 ypos1 pos1 vel1 testTrials

testTrCnt1 = 0; % test trial counter1 for getting through the test dataset
testTrCnt2 = 0; % test trial counter2 for getting through the test phase
for f = 1:uv.fold % increment cross-validation folds
    
    stackBinSCtrain = stackBinSC(~testTrialLogic{f,1},:,:);         % select train trials
    stackBinSCtest  = stackBinSC(testTrialLogic{f,1},:,:);          % select test trials
    PosTrain        = pos1(~testTrialLogic{f,1},uv.behavTimeIdx);   % Pos train trials
    PosTest         = pos1(testTrialLogic{f,1},uv.behavTimeIdx);    % Pos test trials
    
    % Train dataset
    for t=1:size(stackBinSCtrain,1) % increment train trials
        curSpikes  = squeeze(stackBinSCtrain(t,:,:));           % timebin (1ms) x neuron
        timeWinE   = uv.timeStep:uv.timeStep:size(curSpikes,1); % time windows with steps of 20 ms (default)
        curSpikes  = curSpikes(1:timeWinE(end-1),:);            % curtail the last time step
        curPos     = PosTrain(t,:);                             % current trial's position
        curPos     = curPos(:,1:timeWinE(end));
        spikeMtr   = reshape(curSpikes,uv.timeStep,[],size(curSpikes,2)); % get spike count matrices (time x timesteps x neurons)
        train(t,f).spikeCount = squeeze(sum(spikeMtr,1))'; % neuron x timesteps
        pos        = curPos(1,timeWinE);                   % get the position at timesteps
        vel        = (pos(:,2:end)-pos(:,1:end-1))/uv.timeStep; % get the binned velocity
        train(t,f).state = [pos(:,1:end-1); vel];          % get the state matrix (position & velocity)
    end
    clearvars t
    
    % Test dataset
    for t=1:size(stackBinSCtest,1)   % increment test trials
        testTrCnt1 = testTrCnt1 + 1; % count test trials
        curSpikes  = squeeze(stackBinSCtest(t,:,:));            % timebin (1ms) x neuron
        timeWinE   = uv.timeStep:uv.timeStep:size(curSpikes,1); % time windows with steps of 20 ms (default)
        curSpikes  = curSpikes(1:timeWinE(end-1),:);            % curtail the last time step
        curPos     = PosTest(t,:);          % current trial's x-y position
        curPos     = curPos(:,1:timeWinE(end));
        spikeMtr   = reshape(curSpikes,uv.timeStep,[],size(curSpikes,2));    % get spike count matrices (time x timesteps x neurons)
        test(testTrials(testTrCnt1)).spikeCount = squeeze(sum(spikeMtr,1))'; % neuron x timesteps
        pos        = curPos(1,timeWinE);                        % get the position at timesteps
        vel        = (pos(:,2:end)-pos(:,1:end-1))/uv.timeStep; % get the binned velocity
        test(testTrials(testTrCnt1)).state = [pos(:,1:end-1); vel];                 % get the state matrix (position & velocity)
        test(testTrials(testTrCnt1)).statePeaks = max([pos(:,1:end-1); vel],[],2);  
    end
    clearvars t 
    
    %% %%%%%%%%%%%%%%%%%%%% Training Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit A, Pi, V and C
    interStateMtr=0;    % inter-state matrix
    intraStateMtr=0;    % intra-state matrix
    obsStateMtr=0;      % observed state matrix
    obsIntraStateMtr=0; % observed intrastate matrix
    count=0;
    countPool=0;
    
    for t=1:length(train) % increment train trials
        % for parameter A, state model describes how state evolves over time
        stateZ1=train(t).state(:,2:end);     % Zt
        stateZ2=train(t).state(:,1:end-1);   % Zt-1
        tmp1=stateZ1*stateZ2';               % Zt*Zt-1' (4-by-4 matrix)
        interStateMtr=interStateMtr+tmp1;    % Sum Zt*Zt-1'
        tmp2=stateZ2*stateZ2';               % Zt-1*Zt-1' (4-by-4 matrix)
        intraStateMtr=intraStateMtr+tmp2;    % Sum Zt-1*Zt-1'
        % for parameter C, observation model describes how observation relates to the state
        stateZ=train(t).state;       % State (position, velocity)
        obserX=train(t).spikeCount;  % Observation (spike data)
        tmp3=obserX*stateZ';                    % Xt*Zt' (#unit-by-4 matrix)
        obsStateMtr=obsStateMtr+tmp3;           % Sum Xt*Zt'
        tmp4=stateZ*stateZ';                    % Zt*Zt'
        obsIntraStateMtr=obsIntraStateMtr+tmp4; % Sum Zt*Zt'
        % for parameter Pi and V
        count=count+1;
        poolZStart(:,count)=train(t).state(:,1);
    end
    clearvars t tmp* count
    
    A=(interStateMtr)*inv(intraStateMtr);   % Slope for the state model, which defines how state evolves over time
    C=(obsStateMtr)*inv(obsIntraStateMtr);  % Slope for the observation model, which defines how observation relates to the state
    Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
    V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
    
    % Fit Q and R (covariances of the state and observation model)
    sumQMtr=0;
    sumRMtr=0;
    count=0;
    for t=1:length(train) % increment train trials
        % for parameter Q, covariance of the state model
        stateZ1=train(t).state(:,2:end);   % Zt
        stateZ2=train(t).state(:,1:end-1); % Zt-1
        tmp5=stateZ1-A*stateZ2;            % Zt-A*Zt-1
        tmp5=tmp5*tmp5';                   % (Zt-AZt-1)(Zt-AZt-1)'
        count=count+size(stateZ1,2);
        sumQMtr=sumQMtr+tmp5;              % Sum (Zt-AZt-1)(Zt-AZt-1)'
        % for parameter R, covariance of the observation model
        stateZ=train(t).state;             % Zt
        obserX=train(t).spikeCount;        % Xt (Neuron-by-timebins)
        tmp6=obserX-C*stateZ;              % Xt-C*Zt
        tmp6=tmp6*tmp6';                   % (Xt-C*Zt)(Xt-C*Zt)'
        sumRMtr=sumRMtr+tmp6;              % Sum (Xt-C*Zt)(Xt-C*Zt)'
    end
    clearvars t tmp*
    
    Q=sumQMtr/(count);                % count 1/(t-1), count should correspond to 2 to T
    R=sumRMtr/(count+length(train));  % count 1/t, count should correspond to 1 to T
    
    % store fitted parameters of each fold
    kfParams(f).A  = A;
    kfParams(f).Pi = Pi;
    kfParams(f).V  = V;
    kfParams(f).C  = C;
    kfParams(f).Q  = Q;
    kfParams(f).R  = R;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% Test Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t=1:size(stackBinSCtest,1)  % increment test trials
        testTrCnt2 = testTrCnt2 + 1;
        numbBin=size(test(testTrials(testTrCnt2)).state,2); % # of bins
        % Initialization
        mu=Pi;   % sample mean
        sigma=V; % sample covariance
        clearvars curEstStateMean curEstStateCov
        for b=1:numbBin % increment timebins (varies trial-by-trial)
            % One-step Prediction
            mu=A*mu;            % A is the coefficient matrix for the state model that maps the previous states to current states
            sigma=A*sigma*A'+Q; % Q is the covariance of the state model
            % compute the Kalman gain
            K=sigma*C'*inv(C*sigma*C'+R); % Kalman gain
            % Measurement update
            curObservX=test(testTrials(testTrCnt2)).spikeCount(:,b); % take the current spike-count bin
            mu=mu+K*(curObservX-C*mu);
            sigma=sigma-K*C*sigma;
            
            curEstStateMean(:,b) =mu;    % current estimation (bin-by-bin) of the state (mean)
            curEstStateCov(:,:,b)=sigma; % current estimation (bin-by-bin) of the state (covariance)
        end
        clearvars b
        test(testTrials(testTrCnt2)).estStateMean=curEstStateMean;  % estimation of the state trajectory (mean)
        test(testTrials(testTrCnt2)).estStateCov =curEstStateCov;   % estimation of the state trajectory (covariance)
        test(testTrials(testTrCnt2)).estStatePeaks=max(curEstStateMean,[],2); % peak position and velocity 
        
    end   
    clearvars t 
end
clearvars -except binSC stackBinSC uv testTrialLogic xpos1 ypos1 pos1 vel1 testTrials train test kfParams

save('KF_DMSonly_090117','train','test','kfParams')

%% %%%%%%%%%%%%%%%%%%%% Plot the selected trials %%%%%%%%%%%%%%%%%%%%%%%%%%
% selTrialList=[ 1, 2, 3, 4, 5, 6, 7, 8 ];
% 
% for t=1:length(selTrialList) % # of selected trials
% 
%     curTrial=test(selTrialList(t)); % (Trial,Direction)
%     figure(t);
%     
%     % 1-d trajectory pos
%     subplot(1,2,1);
%     plot(curTrial.state(1,:),'-ok','LineWidth',2); % actual state trajectory
%     hold on;
%     plot(curTrial.estStateMean(1,:),'-or','LineWidth',2); % estimated state trajectory
%     errorbar(curTrial.estStateMean(1,:),squeeze(curTrial.estStateCov(1,1,:))','r','LineWidth',2)
%     xlabel('Time'); ylabel('Pos');
%     
%     % 1-d trajectory vel
%     subplot(1,2,2); % horizontal position
%     plot(curTrial.state(2,:),'-ok','LineWidth',2); % actual state trajectory
%     hold on;
%     plot(curTrial.estStateMean(2,:),'-or','LineWidth',2); % estimated state trajectory
%     errorbar(curTrial.estStateMean(2,:),squeeze(curTrial.estStateCov(2,2,:))','r','LineWidth',2)
%     xlabel('Time'); ylabel('Vel');
% 
% end
% clearvars t

%% %%%%%%%% Plot the continuous actual/estimated trajectories %%%%%%%%%%%%%
% frameCnt = 0;     % frame counter
% conState    = []; % continuous state
% conEstState = []; % continuous estimated state
% 
% v = VideoWriter('reachTraj.mp4'); % create a video object
% open(v);
% hold on;
% for t = 1:20 % length(test) % increment trials (reaches)
%     for b = 1:length(test(t).state(1,:))  % increment bins (timebins: e.g. 40 ms)
%         frameCnt     = frameCnt + 1;      % frame count
%         tmpState     = test(t).state(1,b);         % current actual position trajectory
%         tmpEstState  = test(t).estStateMean(1,b);  % current estimated position trajectory
%         conState    = [conState,tmpState];         % append actual state
%         conEstState = [conEstState, tmpEstState];  % append estimated state
%         
%         plot(conState,'-k','LineWidth',2);         % plot the appended actual state
%         plot(conEstState,'-r','LineWidth',2);      % plot the appended estimated state
%         
%         % this part is just to keep the x-axis of the plot constantly at
%         % three times of one reach trajectory. 
%         if length(conState) <= 3*length(test(t).state(1,:))  
%             xlim([0 3*length(test(t).state(1,:))]); %
%         elseif length(conState) > 3*length(test(t).state(1,:)) 
%             xlim([length(conState)-3*length(test(t).state(1,:))  length(conState)]); %
%         end
%         
%         % get the current frame to construct the video
%         M(frameCnt)  = getframe;   % register each frame 
%         writeVideo(v,M(frameCnt))  % write the frame to the video file
%         fprintf('processed %d\n', frameCnt)
%     end 
% end
% hold off; clearvars t b frameCnt
% close(v);

%axes('Position',[0 0 1 1])

% replay the movie
%movie(M,1)

%% metric quantifying match between estimated and actual reach trajectories (amplitude)
% correlation between peaks of estimated and actual reach trajectories
testCell  = squeeze(struct2cell(test)); % convert test structure to a cell
testActPeaks = cell2mat(testCell(3,:)); % get actual position peak values
testEstPeaks = cell2mat(testCell(6,:)); % get estimated position peak values
[posCorrR,posCorrP] = corr(testActPeaks(1,:)', testEstPeaks(1,:)'); % the correlation between actual and estimated max positions 
%scatter(testActPeaks(1,:)',testEstPeaks(1,:)');                     % scatterplot to reveal the correlation between actual and estimated max positions  
[velCorrR,velCorrP] = corr(testActPeaks(2,:)', testEstPeaks(2,:)'); % the correlation between actual and estimated max positions 
%scatter(testActPeaks(2,:)',testEstPeaks(2,:)');                     % scatterplot to reveal the correlation between actual and estimated max positions  




