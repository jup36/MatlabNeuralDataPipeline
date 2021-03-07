function [bTj, bTjStimLaser, bTjPStimLaser, lTj] = stimReachKinematicsInDetailForCollectStmLaserCorrect( filePath, pseudoLaserLogic )
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
load(fullfile('/Volumes/Beefcake/Junchol_Data/oldJoystickCalib/jsCalibData.mat'),'mmConv*'); % voltage to mm conversion coefficients from the calibration data

% conversion from voltage to mm
reach0mm = reach0.*mmConvR0;
xPosmm = positionData(1,:).*mmConvX;
yPosmm = positionData(2,:).*mmConvY;

%% get all the rewarded trials and check the coincidence of laser stim   
ts.stmLaser = [ts.stmLaser(1), ts.stmLaser(ts.stmLaser(2:end)-ts.stmLaser(1:end-1)>4000)];
rewardArrC  = arrayfun(@(x) x-2000:x, ts.reward, 'UniformOutput',false); % pre-reward period
rewardStimI = cellfun(@(y) ismember(ts.stmLaser,y), rewardArrC, 'UniformOutput', false); % index for stim delivery during the pre-reward period 
rewardStimIdx = cellfun(@(y) find(y==true), rewardStimI, 'UniformOutput', false ); 
rewardStim = cell2mat(cellfun(@(y) sum(ismember(ts.stmLaser,y)), rewardArrC, 'UniformOutput', false)); % logical for stim delivery during the pre-reward period 

% rewardPStimI = cellfun(@(y) ismember(ts.pseudoLaser,y), rewardArrC, 'UniformOutput', false); % index for the rewarded trials with pseudo laser stim
% rewardPStimIdx = cellfun(@(y) find(y==true), rewardPStimI, 'UniformOutput', false ); 
% rewardPStim = cell2mat(cellfun(@(y) sum(ismember(ts.pseudoLaser,y)), rewardArrC, 'UniformOutput', false)); % logical for stim delivery during the pre-reward period 
% rewardStim(logical(rewardPStim))=2; % put 2 for the pseudoStim trials

% sum(rewardStim&rewardPStim) % just a sanity check for the stim and pseudostim trials being mutually exclusive 

%% get the event windows
bTS.rwd = ts.reward;
bTS.rwdBinE1ms  = -3000:2000; % 1ms bin
bTS.rwdBinE50ms = -3000:50:2000-1; % 50 ms bin
bTS.rwdReachWin = [-2000 1000];  % reach window relative (aligned) to Reward
bTS.rwdLickWin  = [0 2500];      % reward window relative to Reward
bTS.rwdReachBins = bTS.rwdBinE50ms>=bTS.rwdReachWin(1) & bTS.rwdBinE50ms<=bTS.rwdReachWin(2);
bTS.rwdLickBins  = bTS.rwdBinE50ms>=bTS.rwdLickWin(1) & bTS.rwdBinE50ms<=bTS.rwdLickWin(2);

lickBin = zeros(1,length(lick)); % licks
lickBin(ts.lick) = 1; % licks

%% get reward-aligned behavioral kinematic data
gaussianSigma = 2; 
gaussianKernel = TNC_CreateGaussian(gaussianSigma*15,gaussianSigma,gaussianSigma*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

valNoStimRwdTCnt = 0; 
valRwdTs = zeros(length(ts.reward),1);

%% run it with the data aligned to the reward delivery
xPosmmRwd = arrayfun(@(a) xPosmm(a-3000:a+2000-1), ts.reward(1:end-1),'UniformOutput',false);
yPosmmRwd = arrayfun(@(a) yPosmm(a-3000:a+2000-1), ts.reward(1:end-1),'UniformOutput',false);

% rezero x, y trajectories
xPosmmRwd0 = cellfun(@(a) a-median(a(1:500)), xPosmmRwd,'UniformOutput', false); % re-zeroed x traj
yPosmmRwd0 = cellfun(@(a) a-median(a(1:500)), yPosmmRwd,'UniformOutput', false); % re-zeroed y traj

% get the reach trajectories summed across x and y
xPosmmRwd0Sq = cellfun(@(a) a.^2,xPosmmRwd0,'UniformOutput',false);
yPosmmRwd0Sq = cellfun(@(a) a.^2,yPosmmRwd0,'UniformOutput',false);
reachXYrwd = cellfun(@(a,b) sqrt(a+b), xPosmmRwd0Sq, yPosmmRwd0Sq,'UniformOutput',false);

% get lick times relative to each pseudo laser
rwdrchTimeArrC = arrayfun(@(a) a-1500:a-1000, ts.reward(1:end),'UniformOutput',false);  
lTj.lickRelRwdRchC = cellfun(@(a) ts.lick(ismember(ts.lick,a))-a(1), rwdrchTimeArrC,'UniformOutput',false); 
lTj.lickCountRwdRch = cell2mat(cellfun(@length, lTj.lickRelRwdRchC,'UniformOutput',false)); 

for t = 1:length(reachXYrwd) % increment trials, take the trial-by-trial position/velocity data and bin them to match the neural trjectories (e.g. 50 ms)
    if ts.reward(t)+bTS.rwdBinE1ms(end)<=length(reach0mm) && ts.reward(t)+bTS.rwdBinE1ms(1)>0 && rewardStim(t)~=1
        valRwdTs(t,1) = true;
        valNoStimRwdTCnt = valNoStimRwdTCnt+1; 
        timeWin = ts.reward(t)+bTS.rwdBinE1ms; % the time window, e.g. -3 to 2 sec relative to the behavioral timestamp
        bTj(valNoStimRwdTCnt).reachPos = smooth(binAvg1msSpkCountMat(reachXYrwd{t},50,50),3)'; % get decimated behavioral trjectories on the same timescale of the neural trjectories
        bTj(valNoStimRwdTCnt).reachPosReachWin = bTj(valNoStimRwdTCnt).reachPos(bTS.rwdReachBins);   
        bTj(valNoStimRwdTCnt).reachVel = smooth(diff([bTj(valNoStimRwdTCnt).reachPos(1) bTj(valNoStimRwdTCnt).reachPos]),3)'.*(1000/50); % get reach velocities (from reach0mm)
        bTj(valNoStimRwdTCnt).reachVelReachWin = bTj(valNoStimRwdTCnt).reachVel(bTS.rwdReachBins); % get reach velocities (from reach0mm)
        bTj(valNoStimRwdTCnt).maxReachPos = max(bTj(valNoStimRwdTCnt).reachPosReachWin); % max reach position within the time window of interest
        bTj(valNoStimRwdTCnt).maxReachVel = max(bTj(valNoStimRwdTCnt).reachVelReachWin); % max reach velocity within the time window of interest
        bTj(valNoStimRwdTCnt).trialId = t;
        bTj(valNoStimRwdTCnt).lick = bin1msSpkCountMat(lickBin(timeWin),50, 50);   % get binned lick counts using the bin1msSpkCountMat
        bTj(valNoStimRwdTCnt).lickTrace = conv(bTj(valNoStimRwdTCnt).lick, gaussianKernel, 'same')*(1000); % smoothing with a Gaussian kernel
        bTj(valNoStimRwdTCnt).lickCount = sum(bTj(valNoStimRwdTCnt).lick(bTS.rwdLickBins)); % lick Counts within the lickTimeBins               
    end
end
clearvars t

%% get behavioral kinematic data aligned to laser stims that led to reward
xPosmmStim = arrayfun(@(a) xPosmm(a-2000:a+3000-1), ts.stmLaser(1:end-1),'UniformOutput',false);
yPosmmStim = arrayfun(@(a) yPosmm(a-2000:a+3000-1), ts.stmLaser(1:end-1),'UniformOutput',false);

% rezero x, y trajectories
xPosmmStim0 = cellfun(@(a) a-median(a(1:500)), xPosmmStim,'UniformOutput', false); % re-zeroed x traj
yPosmmStim0 = cellfun(@(a) a-median(a(1:500)), yPosmmStim,'UniformOutput', false); % re-zeroed y traj

% get the reach trajectories summed across x and y
xPosmmStim0Sq = cellfun(@(a) a.^2,xPosmmStim0,'UniformOutput',false);
yPosmmStim0Sq = cellfun(@(a) a.^2,yPosmmStim0,'UniformOutput',false);
reachXYStim = cellfun(@(a,b) sqrt(a+b), xPosmmStim0Sq, yPosmmStim0Sq,'UniformOutput',false);

% get lick times relative to each stim laser
stimTimeArrC = arrayfun(@(a) a:a+500, ts.stmLaser,'UniformOutput',false);  
lTj.lickRelStimC = cellfun(@(a) ts.lick(ismember(ts.lick,a))-a(1), stimTimeArrC,'UniformOutput',false); 
lTj.lickCountStim = cell2mat(cellfun(@length, lTj.lickRelStimC,'UniformOutput',false)); 

bTS.stmLaserRwd = ts.stmLaser; %(cell2mat(rewardStimIdx)); % laser stim trials that led to reward  
bTS.rchBinE1ms  = -2000:3000; % 1ms bin
bTS.rchBinE50ms = -2000:50:3000-1; % 50 ms bin
bTS.rchReachWin = [-1000 2000];  % reach window relative to reachStart
bTS.rchLickWin  = [500 3000];    % reward window relatie to reachStart
bTS.rchReachBins = bTS.rchBinE50ms>=bTS.rchReachWin(1) & bTS.rchBinE50ms<=bTS.rchReachWin(2);
bTS.rchLickBins  = bTS.rchBinE50ms>=bTS.rchLickWin(1) & bTS.rchBinE50ms<=bTS.rchLickWin(2);

valStmLaserRwd = zeros(length(reachXYStim),1); 
for t = 1:length(reachXYStim) 
    if bTS.stmLaserRwd(t)+bTS.rchBinE1ms(end)<=length(reach0mm) && bTS.stmLaserRwd(t)+bTS.rchBinE1ms(1)>0
        valStmLaserRwd(t,1) = 1; 
        timeWin = bTS.stmLaserRwd(t)+bTS.rchBinE1ms; % the time window relative to reachStart, -2 to 3 sec
        bTjStimLaser(t).reachPos = smooth(binAvg1msSpkCountMat(reachXYStim{t},50,50),3)'; % get reach position
        bTjStimLaser(t).reachPosReachWin = bTjStimLaser(t).reachPos(bTS.rchReachBins); % reach Position within the reach window
        bTjStimLaser(t).reachVel = smooth(diff([bTjStimLaser(t).reachPos(1) bTjStimLaser(t).reachPos]),3)'.*(1000); % get reach velocities 
        bTjStimLaser(t).reachVelReachWin = bTjStimLaser(t).reachVel(bTS.rchReachBins); % reach velocity within the reach window
        bTjStimLaser(t).maxReachPos = max(bTjStimLaser(t).reachPosReachWin); % max reach position within the time window of interest
        bTjStimLaser(t).maxReachVel = max(bTjStimLaser(t).reachVelReachWin); % max reach velocity within the time window of interest        
        bTjStimLaser(t).trialId = t;
        bTjStimLaser(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),50, 50);   % get binned lick counts using the bin1msSpkCountMat
        bTjStimLaser(t).lickTrace = conv(bTjStimLaser(t).lick, gaussianKernel, 'same')*(1000); % smoothing with a Gaussian kernel
        bTjStimLaser(t).lickCount = sum(bTjStimLaser(t).lick(bTS.rchLickBins)); % lick Counts within the lickTimeBins               
    end
    
end
clearvars t

%% get behavioral kinematic data aligned to pseudo-laser stims that led to reward
if pseudoLaserLogic
    if contains(filePath,'IT06')
        ts.pseudoLaser = ts.pseudoLaser-1000; 
    end

    valPseudoLaserTs = ts.pseudoLaser(ts.pseudoLaser<length(positionData)-3000 & ts.pseudoLaser>2000);
    
    xPosmmPStim = arrayfun(@(a) xPosmm(a-2000:a+3000-1), valPseudoLaserTs(1:end-1),'UniformOutput',false);
    yPosmmPStim = arrayfun(@(a) yPosmm(a-2000:a+3000-1), valPseudoLaserTs(1:end-1),'UniformOutput',false);
    
    % rezero x, y trajectories
    xPosmmPStim0 = cellfun(@(a) a-median(a(1:500)), xPosmmPStim,'UniformOutput', false); % re-zeroed x traj
    yPosmmPStim0 = cellfun(@(a) a-median(a(1:500)), yPosmmPStim,'UniformOutput', false); % re-zeroed y traj
    
    % get the reach trajectories summed across x and y
    xPosmmPStim0Sq = cellfun(@(a) a.^2,xPosmmPStim0,'UniformOutput',false);
    yPosmmPStim0Sq = cellfun(@(a) a.^2,yPosmmPStim0,'UniformOutput',false);
    reachXYPStim = cellfun(@(a,b) sqrt(a+b), xPosmmPStim0Sq, yPosmmPStim0Sq,'UniformOutput',false);
    
    % get lick times relative to each pseudo laser
    pstimTimeArrC = arrayfun(@(a) a:a+500, ts.pseudoLaser,'UniformOutput',false);  
    lTj.lickRelPstimC = cellfun(@(a) ts.lick(ismember(ts.lick,a))-a(1), pstimTimeArrC,'UniformOutput',false); 
    lTj.lickCountPstim = cell2mat(cellfun(@length, lTj.lickRelPstimC,'UniformOutput',false)); 
    
    bTS.pstmLaserRwd = ts.pseudoLaser; %(cell2mat(rewardPStimIdx)); % pseudo-laser stim trials that led to reward
    valPStmLaserRwd = zeros(length(reachXYPStim ),1);
    for t = 1:length(reachXYPStim)
        if bTS.pstmLaserRwd(t)+bTS.rchBinE1ms(end)<=length(reach0mm) && bTS.pstmLaserRwd(t)+bTS.rchBinE1ms(1)>0
            valPStmLaserRwd(t,1) = 1;
            timeWin = bTS.pstmLaserRwd(t)+bTS.rchBinE1ms; % the time window relative to reachStart, -2 to 3 sec
            bTjPStimLaser(t).reachPos = smooth(binAvg1msSpkCountMat(reachXYPStim{t},50,50),3)'; % get reach position
            bTjPStimLaser(t).reachPosReachWin = bTjPStimLaser(t).reachPos(bTS.rchReachBins); % reach Position within the reach window
            bTjPStimLaser(t).reachVel = smooth(diff([bTjPStimLaser(t).reachPos(1) bTjPStimLaser(t).reachPos]),3)'.*(1000); % get reach velocities
            bTjPStimLaser(t).reachVelReachWin = bTjPStimLaser(t).reachVel(bTS.rchReachBins); % reach velocity within the reach window
            bTjPStimLaser(t).maxReachPos = max(bTjPStimLaser(t).reachPosReachWin); % max reach position within the time window of interest
            bTjPStimLaser(t).maxReachVel = max(bTjPStimLaser(t).reachVelReachWin); % max reach velocity within the time window of interest
            bTjPStimLaser(t).trialId = t;
            bTjPStimLaser(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),50, 50);  % get binned lick counts using the bin1msSpkCountMat
            bTjPStimLaser(t).lickTrace = conv(bTjPStimLaser(t).lick, gaussianKernel, 'same')*(1000); % smoothing with a Gaussian kernel
            bTjPStimLaser(t).lickCount = sum(bTjPStimLaser(t).lick(bTS.rchLickBins)); % lick Counts within the lickTimeBins
        end
    end
    clearvars t
    
else
    bTjPStimLaser = []; 
end

%% plot the reward-aligned behavioral kinematic  data
% c = cellfun(@num2str, num2cell(rewardStim(valRwdTs==1)), 'UniformOutput', false); % class input for color coding in gramm
% c = cellfun(@(y) strrep(y,'0','noStim'), c, 'UniformOutput', false); 
% c = cellfun(@(y) strrep(y,'1','Stim'), c, 'UniformOutput', false); 
% c = cellfun(@(y) strrep(y,'2','p-Stim'), c, 'UniformOutput', false); 
% 
% clear g
% 
% xAxis = bTS.rwdBinE50ms(bTS.rwdReachBins); 
% 
% g(1,1)=gramm('x',xAxis,'y',{bTj(valRwdTs==1).reachPosReachWin},'color',c);
% g(1,1).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,1).set_names('x','Time (ms)','y','Position a.u.');
% g(1,1).set_title('Position');  
% 
% g(1,2)=gramm('x',xAxis,'y',{bTj(valRwdTs==1).reachVelReachWin},'color',c);
% g(1,2).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,2).set_names('x','Time (ms)','y','Velocity a.u.');
% g(1,2).set_title('Velocity');  
% 
% g(1,3)=gramm('x',bTS.rwdBinE50ms,'y',{bTj(valRwdTs==1).lickTrace},'color',c);
% g(1,3).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,3).set_names('x','Time (ms)','y','lick Counts per sec');
% g(1,3).set_title('Lick Count');  
% 
% figHandle = figure('Position',[100 100 800 350]);
% g.set_title('Kinematics aligned to Reward'); 
% g.draw();

%print( fullfile(filePath,'Figure',strcat(saveNameTag,'rewardedStimVsNoStimVsPStimKinematics_alignToReward')), '-dpdf', '-bestfit')

%% plot the stimLaser and p-stimLaser aligned behavioral kinematic data
% c = cellfun(@num2str, num2cell([ones(length(bTjStimLaser(valStmLaserRwd==1)),1); zeros(length(bTjPStimLaser(valPStmLaserRwd==1)),1)]'),'UniformOutput', false); % class input for color coding in gramm
% c = cellfun(@(y) strrep(y,'1','Stim'), c, 'UniformOutput', false); 
% c = cellfun(@(y) strrep(y,'0','p-Stim'), c, 'UniformOutput', false); 
% clear g
% 
% xAxis = bTS.rchBinE50ms(bTS.rchReachBins); 
% 
% g(1,1)=gramm('x',xAxis,'y',[{bTjStimLaser(valStmLaserRwd==1).reachPosReachWin}, {bTjPStimLaser(valPStmLaserRwd==1).reachPosReachWin}],'color',c);
% g(1,1).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,1).set_names('x','Time (ms)','y','Position a.u.');
% g(1,1).set_title('Position'); 
% 
% g(1,2)=gramm('x',xAxis,'y',[{bTjStimLaser(valStmLaserRwd==1).reachVelReachWin}, {bTjPStimLaser(valPStmLaserRwd==1).reachVelReachWin}],'color',c);
% g(1,2).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,2).set_names('x','Time (ms)','y','Velocity a.u.');
% g(1,2).set_title('Velocity'); 
% 
% % g(1,3)=gramm('x',bTS.rchBinE50ms,'y',[{bTjStimLaser(valStmLaserRwd==1).lickTrace}, {bTjPStimLaser(valPStmLaserRwd==1).lickTrace}],'color',c); %, {bTjPStimLaser(valPStmLaserRwd==1).lickTrace}],'color',c);
% % g(1,3).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% % g(1,3).set_names('x','Time (ms)','y','lick Counts per sec');
% % g(1,3).set_title('Lick Count'); 
% 
% figHandle = figure('Position',[100 100 800 350]);
% g.set_title('Kinematics aligned to Laser vs pseudoLaser'); 
% g.draw();
% 
% print( fullfile(filePath,'Figure',strcat(saveNameTag,'stimVSpstimKinematics_alignToStimOn')), '-dpdf', '-bestfit')

%% plot the stimLaser aligned behavioral kinematic data
% clear g
% 
% xAxis = bTS.rchBinE50ms(bTS.rchReachBins); 
% 
% g(1,1)=gramm('x',xAxis,'y',{bTjStimLaser(valStmLaserRwd==1).reachPosReachWin},'color',c);
% g(1,1).geom_line(); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,1).set_names('x','Time (ms)','y','Position a.u.');
% 
% g(1,2)=gramm('x',xAxis,'y',{bTjStimLaser(valStmLaserRwd==1).reachVelReachWin},'color',c);
% g(1,2).geom_line(); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,2).set_names('x','Time (ms)','y','Velocity a.u.');
% 
% g(1,3)=gramm('x',bTS.rchBinE50ms,'y',{bTjStimLaser(valStmLaserRwd==1).lickTrace},'color',c);
% g(1,3).geom_line(); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,3).set_names('x','Time (ms)','y','lick Counts per sec');
% 
% figHandle = figure('Position',[100 100 800 350]);
% g.draw();
% 
% print( fullfile(filePath,'Figure','stimLaserRewardBehavioralKinematicsIndividual'), '-dpdf', '-bestfit')

%% save Data
% save(fullfile(filePath,strcat(saveNameTag,'_stimReachKinematicsInDetail')), 'bTj', 'bTS', 'bTjStimLaser')

end