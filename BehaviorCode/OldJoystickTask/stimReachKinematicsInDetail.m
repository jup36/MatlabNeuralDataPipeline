function [bTj, bTjStimLaser, bTjPStimLaser] = stimReachKinematicsInDetail( filePath, saveNameTag )
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
load(fullfile(filePath,'BehVariables.mat'),'ts','reach0','lick')

%% get all the rewarded trials and check the coincidence of laser stim   
rewardArrC  = arrayfun(@(x) x-2000:x, ts.reward, 'UniformOutput',false); % pre-reward period
rewardStimI = cellfun(@(y) ismember(ts.stmLaser,y), rewardArrC, 'UniformOutput', false); % index for stim delivery during the pre-reward period 
rewardStimIdx = cellfun(@(y) find(y==true), rewardStimI, 'UniformOutput', false ); 
rewardStim = cell2mat(cellfun(@(y) sum(ismember(ts.stmLaser,y)), rewardArrC, 'UniformOutput', false)); % logical for stim delivery during the pre-reward period 

rewardPStimI = cellfun(@(y) ismember(ts.pseudoLaser,y), rewardArrC, 'UniformOutput', false); % index for the rewarded trials with pseudo laser stim
rewardPStimIdx = cellfun(@(y) find(y==true), rewardPStimI, 'UniformOutput', false ); 
rewardPStim = cell2mat(cellfun(@(y) sum(ismember(ts.pseudoLaser,y)), rewardArrC, 'UniformOutput', false)); % logical for stim delivery during the pre-reward period 
rewardStim(logical(rewardPStim))=2; % put 2 for the pseudoStim trials

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

valRwdTs = zeros(length(ts.reward),1); 
for t = 1:length(ts.reward) % increment trials, take the trial-by-trial position/velocity data and bin them to match the neural trjectories (e.g. 50 ms)
    if ts.reward(t)+bTS.rwdBinE1ms(end)<=length(reach0) && ts.reward(t)+bTS.rwdBinE1ms(1)>0
        valRwdTs(t,1) = true; 
        timeWin = ts.reward(t)+bTS.rwdBinE1ms; % the time window, e.g. -3 to 2 sec relative to the behavioral timestamp
        bTj(t).reachPos = smooth(binAvg1msSpkCountMat(reach0(timeWin),50,50),3)'; % get decimated behavioral trjectories on the same timescale of the neural trjectories
        bTj(t).reachPosReachWin = bTj(t).reachPos(bTS.rwdReachBins);   
        bTj(t).reachVel = smooth(diff([bTj(t).reachPos(1) bTj(t).reachPos]),3)'.*1000; % get reach velocities (from reach0)
        bTj(t).reachVelReachWin = bTj(t).reachVel(bTS.rwdReachBins); % get reach velocities (from reach0)
        bTj(t).maxReachPos = max(bTj(t).reachPosReachWin); % max reach position within the time window of interest
        bTj(t).maxReachVel = max(bTj(t).reachVelReachWin); % max reach velocity within the time window of interest
        bTj(t).trialId = t;
        bTj(t).lick = bin1msSpkCountMat(lickBin(timeWin),50, 50);   % get binned lick counts using the bin1msSpkCountMat
        bTj(t).lickTrace = conv(bTj(t).lick, gaussianKernel, 'same')*(1000/50); % smoothing with a Gaussian kernel
        bTj(t).lickCount = sum(bTj(t).lick(bTS.rwdLickBins)); % lick Counts within the lickTimeBins               
    end
end
clearvars t

%% get behavioral kinematic data aligned to laser stims that led to reward
bTS.stmLaserRwd = ts.stmLaser; %(cell2mat(rewardStimIdx)); % laser stim trials that led to reward  
bTS.rchBinE1ms  = -2000:3000; % 1ms bin
bTS.rchBinE50ms = -2000:50:3000-1; % 50 ms bin
bTS.rchReachWin = [-1000 2000];  % reach window relative to reachStart
bTS.rchLickWin  = [500 3000];    % reward window relatie to reachStart
bTS.rchReachBins = bTS.rchBinE50ms>=bTS.rchReachWin(1) & bTS.rchBinE50ms<=bTS.rchReachWin(2);
bTS.rchLickBins  = bTS.rchBinE50ms>=bTS.rchLickWin(1) & bTS.rchBinE50ms<=bTS.rchLickWin(2);

valStmLaserRwd = zeros(length(bTS.stmLaserRwd),1); 
for t = 1:length(bTS.stmLaserRwd) 
    if bTS.stmLaserRwd(t)+bTS.rchBinE1ms(end)<=length(reach0) && bTS.stmLaserRwd(t)+bTS.rchBinE1ms(1)>0
        valStmLaserRwd(t,1) = 1; 
        timeWin = bTS.stmLaserRwd(t)+bTS.rchBinE1ms; % the time window relative to reachStart, -2 to 3 sec
        bTjStimLaser(t).reachPos = smooth(binAvg1msSpkCountMat(reach0(timeWin),50,50),3)'; % get reach position
        bTjStimLaser(t).reachPosReachWin = bTjStimLaser(t).reachPos(bTS.rchReachBins); % reach Position within the reach window
        bTjStimLaser(t).reachVel = smooth(diff([bTjStimLaser(t).reachPos(1) bTjStimLaser(t).reachPos]),3)'.*1000; % get reach velocities 
        bTjStimLaser(t).reachVelReachWin = bTjStimLaser(t).reachVel(bTS.rchReachBins); % reach velocity within the reach window
        bTjStimLaser(t).maxReachPos = max(bTjStimLaser(t).reachPosReachWin); % max reach position within the time window of interest
        bTjStimLaser(t).maxReachVel = max(bTjStimLaser(t).reachVelReachWin); % max reach velocity within the time window of interest        
        bTjStimLaser(t).trialId = t;
        bTjStimLaser(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),50, 50);   % get binned lick counts using the bin1msSpkCountMat
        bTjStimLaser(t).lickTrace = conv(bTjStimLaser(t).lick, gaussianKernel, 'same')*(1000/50); % smoothing with a Gaussian kernel
        bTjStimLaser(t).lickCount = sum(bTjStimLaser(t).lick(bTS.rchLickBins)); % lick Counts within the lickTimeBins               
    end
    
end
clearvars t

%% get behavioral kinematic data aligned to pseudo-laser stims that led to reward
bTS.pstmLaserRwd = ts.pseudoLaser; %(cell2mat(rewardPStimIdx)); % pseudo-laser stim trials that led to reward
valPStmLaserRwd = zeros(length(bTS.pstmLaserRwd),1); 
for t = 1:length(bTS.pstmLaserRwd) 
    if bTS.pstmLaserRwd(t)+bTS.rchBinE1ms(end)<=length(reach0) && bTS.pstmLaserRwd(t)+bTS.rchBinE1ms(1)>0
        valPStmLaserRwd(t,1) = 1; 
        timeWin = bTS.pstmLaserRwd(t)+bTS.rchBinE1ms; % the time window relative to reachStart, -2 to 3 sec
        bTjPStimLaser(t).reachPos = smooth(binAvg1msSpkCountMat(reach0(timeWin),50,50),3)'; % get reach position
        bTjPStimLaser(t).reachPosReachWin = bTjPStimLaser(t).reachPos(bTS.rchReachBins); % reach Position within the reach window
        bTjPStimLaser(t).reachVel = smooth(diff([bTjPStimLaser(t).reachPos(1) bTjPStimLaser(t).reachPos]),3)'.*1000; % get reach velocities 
        bTjPStimLaser(t).reachVelReachWin = bTjPStimLaser(t).reachVel(bTS.rchReachBins); % reach velocity within the reach window
        bTjPStimLaser(t).maxReachPos = max(bTjPStimLaser(t).reachPosReachWin); % max reach position within the time window of interest
        bTjPStimLaser(t).maxReachVel = max(bTjPStimLaser(t).reachVelReachWin); % max reach velocity within the time window of interest        
        bTjPStimLaser(t).trialId = t;
        bTjPStimLaser(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),50, 50);  % get binned lick counts using the bin1msSpkCountMat
        bTjPStimLaser(t).lickTrace = conv(bTjPStimLaser(t).lick, gaussianKernel, 'same')*(1000/50); % smoothing with a Gaussian kernel
        bTjPStimLaser(t).lickCount = sum(bTjPStimLaser(t).lick(bTS.rchLickBins)); % lick Counts within the lickTimeBins               
    end
end
clearvars t

%% plot the reward-aligned behavioral kinematic  data
c = cellfun(@num2str, num2cell(rewardStim(valRwdTs==1)), 'UniformOutput', false); % class input for color coding in gramm
c = cellfun(@(y) strrep(y,'0','noStim'), c, 'UniformOutput', false); 
c = cellfun(@(y) strrep(y,'1','Stim'), c, 'UniformOutput', false); 
c = cellfun(@(y) strrep(y,'2','p-Stim'), c, 'UniformOutput', false); 

clear g

xAxis = bTS.rwdBinE50ms(bTS.rwdReachBins); 

g(1,1)=gramm('x',xAxis,'y',{bTj(valRwdTs==1).reachPosReachWin},'color',c);
g(1,1).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,1).set_names('x','Time (ms)','y','Position a.u.');
g(1,1).set_title('Position');  

g(1,2)=gramm('x',xAxis,'y',{bTj(valRwdTs==1).reachVelReachWin},'color',c);
g(1,2).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,2).set_names('x','Time (ms)','y','Velocity a.u.');
g(1,2).set_title('Velocity');  

g(1,3)=gramm('x',bTS.rwdBinE50ms,'y',{bTj(valRwdTs==1).lickTrace},'color',c);
g(1,3).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,3).set_names('x','Time (ms)','y','lick Counts per sec');
g(1,3).set_title('Lick Count');  

figHandle = figure('Position',[100 100 800 350]);
g.set_title('Kinematics aligned to Reward'); 
g.draw();

print( fullfile(filePath,'Figure',strcat(saveNameTag,'rewardedStimVsNoStimVsPStimKinematics_alignToReward')), '-dpdf', '-bestfit')

%% plot the stimLaser and p-stimLaser aligned behavioral kinematic data
c = cellfun(@num2str, num2cell([ones(length(bTjStimLaser(valStmLaserRwd==1)),1); zeros(length(bTjPStimLaser(valPStmLaserRwd==1)),1)]'),'UniformOutput', false); % class input for color coding in gramm
c = cellfun(@(y) strrep(y,'1','Stim'), c, 'UniformOutput', false); 
c = cellfun(@(y) strrep(y,'0','p-Stim'), c, 'UniformOutput', false); 
clear g

xAxis = bTS.rchBinE50ms(bTS.rchReachBins); 

g(1,1)=gramm('x',xAxis,'y',[{bTjStimLaser(valStmLaserRwd==1).reachPosReachWin}, {bTjPStimLaser(valPStmLaserRwd==1).reachPosReachWin}],'color',c);
g(1,1).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,1).set_names('x','Time (ms)','y','Position a.u.');
g(1,1).set_title('Position'); 

g(1,2)=gramm('x',xAxis,'y',[{bTjStimLaser(valStmLaserRwd==1).reachVelReachWin}, {bTjPStimLaser(valPStmLaserRwd==1).reachVelReachWin}],'color',c);
g(1,2).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,2).set_names('x','Time (ms)','y','Velocity a.u.');
g(1,2).set_title('Velocity'); 

% g(1,3)=gramm('x',bTS.rchBinE50ms,'y',[{bTjStimLaser(valStmLaserRwd==1).lickTrace}, {bTjPStimLaser(valPStmLaserRwd==1).lickTrace}],'color',c); %, {bTjPStimLaser(valPStmLaserRwd==1).lickTrace}],'color',c);
% g(1,3).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
% g(1,3).set_names('x','Time (ms)','y','lick Counts per sec');
% g(1,3).set_title('Lick Count'); 

figHandle = figure('Position',[100 100 800 350]);
g.set_title('Kinematics aligned to Laser vs pseudoLaser'); 
g.draw();

print( fullfile(filePath,'Figure',strcat(saveNameTag,'stimVSpstimKinematics_alignToStimOn')), '-dpdf', '-bestfit')

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
save(fullfile(filePath,strcat(saveNameTag,'_stimReachKinematicsInDetail')), 'bTj', 'bTS', 'bTjStimLaser')

end