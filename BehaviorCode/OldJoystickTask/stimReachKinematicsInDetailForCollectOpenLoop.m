function [bTjStimLaser, bTjPStimLaser, thresholdCrossed] = stimReachKinematicsInDetailForCollectOpenLoop( filePath )
%This function is to compare the behavioral kinematic traces of rewarded
% reaches with VS without the laser perturbation. One critical concern on
% comparing unrewarded reaches was whether the animal was really reaching or
% just made a small deflection that triggered stimulation for whatever reason.
% One way to get around this is just focusing on rewarded reaches where we
% know surely that animals were reaching.

%% load unitTimeTrial mat from 'binSpkCount*.mat' and preprocess it
%S = load(fullfile(filePath, binSpkCntFileName), 'reward');
%S = S.('reward');

%% get behavioral data 'BehVariables.mat'
load(fullfile(filePath,'BehVariables.mat'),'ts','reach0','positionData','lick')
load(fullfile('/Volumes/RAID2/parkj/oldJoystickCalib','jsCalibData'),'mmConv*'); % voltage to mm conversion coefficients from the calibration data
reach0mm = reach0.*mmConvR0;

% conversion from voltage to mm
reach0mm = reach0.*mmConvR0;
xPosmm = positionData(1,:).*mmConvX;
yPosmm = positionData(2,:).*mmConvY;

%% get all the previous reward + 2000ms points  and check the coincidence of laser stim
postRwd2s = ts.reward+2000; % stimOn if it's delivered
postRwdArrC  = arrayfun(@(x) x-1000:x+1000, postRwd2s, 'UniformOutput',false); % post-reward period
postRwdStimI = cellfun(@(y) ismember(ts.stmLaser,y), postRwdArrC, 'UniformOutput', false); % index for stim delivery during the post-reward period
postRwdStimIdx = cellfun(@(y) find(y==true), postRwdStimI, 'UniformOutput', false );
postRwdStim = cell2mat(cellfun(@(y) sum(ismember(ts.stmLaser,y)), postRwdArrC, 'UniformOutput', false)); % logical for stim delivery during the post-reward period

% get the thresholdCrossed (check 1-ms bin-by-bin if reach was initiated by each bin)
% stimPeriodArrC = arrayfun(@(x) x:x+5000, postRwd2s(1:end-1), 'UniformOutput',false); % post-reward period
% stimPeriodRwdI = cellfun(@(y) ismember(ts.reward-400,y), stimPeriodArrC, 'UniformOutput', false);
% for i = 1:length(stimPeriodArrC)
%     thresholdCrossedC{i,1} = zeros(1,length(stimPeriodArrC{1,i}))==1;  
%     if sum(stimPeriodRwdI{i})==1
%         thresholdCrossedC{i,1} = stimPeriodArrC{1,i}>ts.reward(stimPeriodRwdI{i})-400; 
%     end
% end
% clearvars i

stimPeriodArrC = arrayfun(@(x) x:x+5000, postRwd2s(1:end-1), 'UniformOutput',false); % post-reward period
stimPeriodRwdI = cellfun(@(y) ismember(ts.reachStart,y), stimPeriodArrC, 'UniformOutput', false);
for i = 1:length(stimPeriodArrC)
    thresholdCrossedC{i,1} = zeros(1,length(stimPeriodArrC{1,i}))==1;  
    if sum(stimPeriodRwdI{i})==1
        thresholdCrossedC{i,1} = stimPeriodArrC{1,i}>ts.reachStart(stimPeriodRwdI{i}); 
    end
end
clearvars i

thresholdCrossed.bTjStimLaserFull = cell2mat(thresholdCrossedC(logical(postRwdStim(1:end-1)))); 
thresholdCrossed.bTjStimLaser = sum(thresholdCrossed.bTjStimLaserFull,1)./size(thresholdCrossed.bTjStimLaserFull,1); 

thresholdCrossed.bTjPstimLaserFull = cell2mat(thresholdCrossedC(~logical(postRwdStim(1:end-1)))); 
thresholdCrossed.bTjPstimLaser = sum(thresholdCrossed.bTjPstimLaserFull,1)./size(thresholdCrossed.bTjPstimLaserFull,1); 

%% get reward-aligned behavioral kinematic data
gaussianSigma = 1;
gaussianKernel = TNC_CreateGaussian(gaussianSigma*15,gaussianSigma,gaussianSigma*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

lickBin = zeros(1,length(lick)); % licks
lickBin(ts.lick) = 1; % licks

%% get behavioral kinematic data aligned to laser stims
xPosmmStim = arrayfun(@(a) xPosmm(a-1000:a+5000-1), ts.stmLaser(1:end-1),'UniformOutput',false);
yPosmmStim = arrayfun(@(a) yPosmm(a-1000:a+5000-1), ts.stmLaser(1:end-1),'UniformOutput',false);

% rezero x, y trajectories
xPosmmStim0 = cellfun(@(a) a-mean(a(1001:1500)), xPosmmStim,'UniformOutput', false); % re-zeroed x traj
yPosmmStim0 = cellfun(@(a) a-mean(a(1001:1500)), yPosmmStim,'UniformOutput', false); % re-zeroed y traj

% get the reach trajectories summed across x and y
xPosmmStim0Sq = cellfun(@(a) a.^2,xPosmmStim0,'UniformOutput',false);
yPosmmStim0Sq = cellfun(@(a) a.^2,yPosmmStim0,'UniformOutput',false);
reachXYStim = cellfun(@(a,b) sqrt(a+b), xPosmmStim0Sq, yPosmmStim0Sq,'UniformOutput',false);

bTS.stmLaser = ts.stmLaser; %(cell2mat(rewardStimIdx)); % laser stim trials that led to reward
bTS.rchBinE1ms  = 0:5000; % 1ms bin
bTS.rchBinE50ms = 0:50:5000-1; % 50 ms bin

valStmLaserRwd = zeros(length(reachXYStim),1);
for t = 1:length(reachXYStim)
    if ts.stmLaser(t)+bTS.rchBinE1ms(end)<=length(reach0mm) && ts.stmLaser(t)+bTS.rchBinE1ms(1)>0
        valStmLaserRwd(t,1) = 1;
        timeWin = ts.stmLaser(t)+bTS.rchBinE1ms; % the time window relative to reachStart, -2 to 3 sec
        bTjStimLaser(t).reachPos = binAvg1msSpkCountMat(reachXYStim{t},50,50); %smooth(binAvg1msSpkCountMat(reachXYStim{t},50,50),3)'; % get reach position
        bTjStimLaser(t).reachVel = diff([bTjStimLaser(t).reachPos(1) bTjStimLaser(t).reachPos]).*(1000/50); %smooth(diff([bTjStimLaser(t).reachPos(1) bTjStimLaser(t).reachPos]),3)'.*(1000/50); % get reach velocities
        bTjStimLaser(t).maxReachPos = max(bTjStimLaser(t).reachPos); % max reach position within the time window of interest
        bTjStimLaser(t).maxReachVel = max(bTjStimLaser(t).reachVel); % max reach velocity within the time window of interest
        bTjStimLaser(t).trialId = t;
        bTjStimLaser(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),50, 50);   % get binned lick counts using the bin1msSpkCountMat
        bTjStimLaser(t).lickTrace = conv(bTjStimLaser(t).lick, gaussianKernel, 'same')*(1000/50); % smoothing with a Gaussian kernel
        bTjStimLaser(t).lickCount = sum(bTjStimLaser(t).lick(:)); % lick Counts within the lickTimeBins
    end
end
clearvars t

%plot(mean(cell2mat({bTjStimLaser(:).reachPos}'),1))

%% get behavioral kinematic data aligned to pseudo-laser stims that led to reward
ts.pseudoLaser = postRwd2s(~postRwdStim); 
valPseudoLaserTs = ts.pseudoLaser(ts.pseudoLaser<length(positionData)-5000);

xPosmmPStim = arrayfun(@(a) xPosmm(a-1000:a+5000-1), valPseudoLaserTs(1:end-1),'UniformOutput',false);
yPosmmPStim = arrayfun(@(a) yPosmm(a-1000:a+5000-1), valPseudoLaserTs(1:end-1),'UniformOutput',false);

% get the reach trajectories summed across x and y
% rezero x, y trajectories
xPosmmPStim0 = cellfun(@(a) a-mean(a(1:1500)), xPosmmPStim,'UniformOutput', false); % re-zeroed x traj
yPosmmPStim0 = cellfun(@(a) a-mean(a(1:1500)), yPosmmPStim,'UniformOutput', false); % re-zeroed y traj

% get the reach trajectories summed across x and y
xPosmmPStim0Sq = cellfun(@(a) a.^2,xPosmmPStim0,'UniformOutput',false);
yPosmmPStim0Sq = cellfun(@(a) a.^2,yPosmmPStim0,'UniformOutput',false);
reachXYPStim = cellfun(@(a,b) sqrt(a+b), xPosmmPStim0Sq, yPosmmPStim0Sq,'UniformOutput',false);

bTS.pstmLaser = ts.pseudoLaser; %(cell2mat(rewardPStimIdx)); % pseudo-laser stim trials that led to reward
valPStmLaserRwd = zeros(length(reachXYPStim),1);
for t = 1:length(reachXYPStim)
    if ts.pseudoLaser(t)+bTS.rchBinE1ms(end)<=length(reach0mm) && ts.pseudoLaser(t)+bTS.rchBinE1ms(1)>0
        valPStmLaserRwd(t,1) = 1;
        timeWin = ts.pseudoLaser(t)+bTS.rchBinE1ms; % the time window relative to reachStart, -2 to 3 sec
        bTjPStimLaser(t).reachPos = binAvg1msSpkCountMat(reachXYPStim{t},50,50); %smooth(binAvg1msSpkCountMat(reachXYPStim{t},50,50),3)'; % get reach position
        bTjPStimLaser(t).reachVel = diff([bTjPStimLaser(t).reachPos(1) bTjPStimLaser(t).reachPos]).*(1000/50); %smooth(diff([bTjPStimLaser(t).reachPos(1) bTjPStimLaser(t).reachPos]),3)'.*(1000/50); % get reach velocities
        bTjPStimLaser(t).maxReachPos = max(bTjPStimLaser(t).reachPos); % max reach position within the time window of interest
        bTjPStimLaser(t).maxReachVel = max(bTjPStimLaser(t).reachVel); % max reach velocity within the time window of interest
        bTjPStimLaser(t).trialId = t;
        bTjPStimLaser(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),50, 50);  % get binned lick counts using the bin1msSpkCountMat
        bTjPStimLaser(t).lickTrace = conv(bTjPStimLaser(t).lick, gaussianKernel, 'same')*(1000/50); % smoothing with a Gaussian kernel
        bTjPStimLaser(t).lickCount = sum(bTjPStimLaser(t).lick(:)); % lick Counts within the lickTimeBins
    end
end
clearvars t

end