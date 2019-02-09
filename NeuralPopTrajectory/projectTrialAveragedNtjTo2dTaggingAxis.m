function [nTj, bTj] = projectTrialAveragedNtjTo2dTaggingAxis( filePath, tagPCAFileName, binSpkCntFileName, evtName, saveNameTag, varargin )
%This function projects neural activity sorted by behavioral kinematic variables onto  
% the PC axes extracted from tagging responses. Units that passed the criteria for tagging activity PCA
% (pcaResultTag.unitIdx) are used only for projection to get population trajectories within the 2-d tagging state space.  

p = parse_input_projectTrialAveragedNtjTo2dTaggingAxis( filePath, tagPCAFileName, binSpkCntFileName, evtName, saveNameTag, varargin );

%% load tagging activity PCA result
pcaResultTag = load(tagPCAFileName,'pcaResult');
pcaResultTag = pcaResultTag.('pcaResult');
unitIdx = pcaResultTag.unitIdx; % use the pcaResultTag.unitIdx to match the dimension of the tagging activity PCs (to project onto them)
clearvars pcaResult

tagTrjMat  = reshape([pcaResultTag.kern.seqTrain.xpost], pcaResultTag.p.Results.pcaDim, [], length(pcaResultTag.kern.seqTrain)); % neural trajectories in dim x timeBins x trials (e.g. 5x100x114)
avgTagTrjMat = nanmean(tagTrjMat,3); % tag Trjs averaged across whole trials
tagRelativeTimeBins = pcaResultTag.p.Results.timeRange(1):pcaResultTag.p.Results.pcaBinSize:pcaResultTag.p.Results.timeRange(2)-pcaResultTag.p.Results.pcaBinSize; % get timeBins relative to time 0 (e.g. reward) of the current PSTH
tagWin = [0 1000]-pcaResultTag.p.Results.pcaBinSize; % tag window

%% load unitTimeTrial mat from 'binSpkCount*.mat' and preprocess it
S = load(fullfile(filePath, binSpkCntFileName), evtName);
S = S.(evtName);

unitTimeTrialB = bin1msSpkCountMat( S.unitTimeTrial(unitIdx,:,:), p.Results.binSize, p.Results.binSize ); % bin unitTimeTrial mat
% get the mean and std across all 50-ms time bins and trials, then z-score normalize
rsUnitTimeTrial = reshape(unitTimeTrialB, [size(unitTimeTrialB,1), size(unitTimeTrialB,2)*size(unitTimeTrialB,3)]); % reshape
meanPerUnit = nanmean(rsUnitTimeTrial,2); % mean across timeBins and Trials
stdPerUnit = nanstd(rsUnitTimeTrial,0,2); % std across timeBins and Trials
rsUnitTimeTrialBZ=(rsUnitTimeTrial-repmat(meanPerUnit,[1 size(rsUnitTimeTrial,2)]))./repmat(stdPerUnit,[1 size(rsUnitTimeTrial,2)]); % z-score normalize
unitTimeTrialBZ = reshape(rsUnitTimeTrialBZ, [size(unitTimeTrialB,1), size(unitTimeTrialB,2), size(unitTimeTrialB,3)]); % reshape back

%% load unitTimeTrial mat of the stmReach
StmRch = load(fullfile(filePath, binSpkCntFileName), 'stmReach');
StmRch = StmRch.('stmReach'); 

stmRunitTimeTrialB = bin1msSpkCountMat( StmRch.unitTimeTrial(unitIdx,:,:), p.Results.binSize, p.Results.binSize ); % bin unitTimeTrial mat
% get the mean and std across all 50-ms time bins and trials, then z-score normalize
stmRrsUnitTimeTrial = reshape(stmRunitTimeTrialB, [size(stmRunitTimeTrialB,1), size(stmRunitTimeTrialB,2)*size(stmRunitTimeTrialB,3)]); % reshape
stmRmeanPerUnit = nanmean(stmRrsUnitTimeTrial,2); % mean across timeBins and Trials
stmRstdPerUnit = nanstd(stmRrsUnitTimeTrial,0,2); % std across timeBins and Trials
stmRrsUnitTimeTrialBZ=(stmRrsUnitTimeTrial-repmat(stmRmeanPerUnit,[1 size(stmRrsUnitTimeTrial,2)]))./repmat(stmRstdPerUnit,[1 size(stmRrsUnitTimeTrial,2)]); % z-score normalize
stmRunitTimeTrialBZ = reshape(stmRrsUnitTimeTrialBZ, [size(stmRunitTimeTrialB,1), size(stmRunitTimeTrialB,2), size(stmRunitTimeTrialB,3)]); % reshape back
avgStmRchunitTimeTrialBZ = nanmean(stmRunitTimeTrialBZ,3); 

%% get behavioral data 'BehVariables.mat'
load(fullfile(filePath,'BehVariables.mat'),'ts','reach0','lick')
rewardArrC = arrayfun(@(x) x-3000:x, ts.reward, 'UniformOutput',false); % pre-reward period
rewardStim = cell2mat(cellfun(@(y) sum(ismember(ts.stmLaser,y)), rewardArrC, 'UniformOutput', false)); % logical for stim delivery during the pre-reward period 
avgStmRwdUnitTimeTrialBZ = nanmean(unitTimeTrialBZ(:,:,logical(rewardStim)),3); % take the average across trials

%% get the event windows
switch evtName
    case 'reward'
        behTS = ts.reward;
        evtMarkersRelative = [-1000 0]-p.Results.binSize; % time points to mark on the neural trajectories
        reachWin = [-2000 0]; % reach window relative to Reward
        lickWin = [0 2500]; % reward window relative to Reward
    case 'reach'
        behTS = ts.reachStart;
        evtMarkersRelative = [-50 1450]-p.Results.binSize; % reachStart, reward delivery (roughly) to be marked on the neural trajectories
        reachWin = [-500 1000]; % reach window relative to ReachStart
        lickWin = [500 3000]; % reward window relatie to ReachStart
    case 'stmLaser'
        behTS = ts.stmLaser;
        evtMarkersRelative = [-50 1450]-p.Results.binSize; % reachStart, reward delivery (roughly)
        reachWin = [-500 1000]; % reach window relative to ReachStart
        lickWin = [500 3000]; % reward window relatie to ReachStart
    case 'stmReach'
        behTS = ts.stmReachStart;
        evtMarkersRelative = [-50 1450]-p.Results.binSize; % reachStart, reward delivery (roughly)
        reachWin = [-500 1000]; % reach window relative to ReachStart
        lickWin = [500 3000]; % reward window relatie to ReachStart
end

% preprocess the behavioral data
lickBin = zeros(1,length(lick)); % licks
lickBin(ts.lick) = 1; % licks

if length(behTS)~=size(unitTimeTrialBZ,3)
    error('The # of trials in behavioral and neural data do NOT match!')
end

valTrialIdx = zeros(length(behTS),1);

binEdges = S.params.binEdges1ms(1):p.Results.binSize:S.params.binEdges1ms(end)-1;
reachTimeBins = binEdges>=reachWin(1) & binEdges<=reachWin(2);
lickTimeBins = binEdges>=lickWin(1) & binEdges<=lickWin(2);

for t = 1:length(behTS) % increment trials, take the trial-by-trial position/velocity data and bin them to match the neural trjectories (e.g. 50 ms)
    if behTS(t)+S.params.binEdges1ms(end)<=length(reach0) && behTS(t)+S.params.binEdges1ms(1)>0
        valTrialIdx(t)=1;
        timeWin = behTS(t)+S.params.binEdges1ms; % the time window, e.g. -3 to 2 sec relative to the behavioral timestamp
        bTj(t).reachPos = smooth(binAvg1msSpkCountMat(reach0(timeWin),p.Results.binSize,p.Results.binSize),3)'; % get decimated behavioral trjectories on the same timescale of the neural trjectories
        bTj(t).reachVel = smooth(diff([bTj(t).reachPos(1) bTj(t).reachPos]),3)'.*1000; % get reach velocities (from reach0)
        bTj(t).maxReachPos = max(bTj(t).reachPos(reachTimeBins)); % max reach position within the time window of interest
        bTj(t).maxReachVel = max(bTj(t).reachVel(reachTimeBins)); % max reach velocity within the time window of interest
        bTj(t).trialId = t;
        bTj(t).lick = bin1msSpkCountMat(lickBin(timeWin-1),p.Results.binSize, p.Results.binSize); % get binned lick counts using the bin1msSpkCountMat
        bTj(t).lickCount =  sum(bTj(t).lick(lickTimeBins)); % lick Counts within the lickTimeBins
    else
    end
end
clearvars t

% get regressors (position, velocity, lick count)
maxP = [[bTj.maxReachPos];[bTj.trialId]]'; % trial-by-trial max reach positions
maxV = [[bTj.maxReachVel];[bTj.trialId]]'; % trial-by-trial max reach velocities
lickCnt = [[bTj.lickCount];[bTj.trialId]]'; % trial-by-trial lick counts

maxP(:,3) = discretize(maxP(:,1),linspace(min(maxP(:,1)),max(maxP(:,1)),p.Results.numbFolds+1)); % discretize trial-by-trial positions
maxV(:,3) = discretize(maxV(:,1),linspace(min(maxV(:,1)),max(maxV(:,1)),p.Results.numbFolds+1)); % discretize trial-by-trial velocities
lickCnt(:,3) = discretize(lickCnt(:,1),linspace(min(lickCnt(:,1)),max(lickCnt(:,1)),p.Results.numbFolds+1)); % discretize trial-by-trial lick counts

%trFolds = floor(length(bTj)/p.Results.numbFolds); % trial p.Results.numbFolds

%sortMaxPos = sortrows(maxP,-1); % sorted max reach positions in a descending order
%sortMaxVel = sortrows(maxV,-1); % sorted max reach velocities
%sortLickCnt = sortrows(lickCnt,-1);  % sorted trial lick counts

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 1:p.Results.numbFolds
        nTj.sortedFoldMaxPos{f}  = unitTimeTrialBZ(:,:,maxP(:,3)==f); % take the trials of the current fold sorted by the max position
        nTj.sortedFoldMaxVel{f}  = unitTimeTrialBZ(:,:,maxV(:,3)==f); % take the trials of the current fold sorted by the max velocity
        nTj.sortedFoldLickCnt{f} = unitTimeTrialBZ(:,:,lickCnt(:,3)==f); % take the trials of the current fold sorted by the lick counts
end
clearvars f

%% project to tagging activity PCs
projVecs = pcaResultTag.kern.pcDirs(:,p.Results.PCs); % vectors to project neural activity onto
gaussianSigma   = pcaResultTag.kernSDList(1)/p.Results.binSize;  % gaussian std (50ms, as the data are already binned with 50-ms window)
gaussianKernel  = TNC_CreateGaussian(gaussianSigma*15,gaussianSigma,gaussianSigma*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

nTj.projToTagAxesSortByMaxPos = projectNtjToAxesSmooth( nTj.sortedFoldMaxPos, projVecs, gaussianKernel ); 
nTj.projToTagAxesSortByMaxVel = projectNtjToAxesSmooth( nTj.sortedFoldMaxVel, projVecs, gaussianKernel ); 
nTj.projToTagAxesSortByLickCnt = projectNtjToAxesSmooth( nTj.sortedFoldLickCnt, projVecs, gaussianKernel ); 

nTj.projToTagAxesEvtAvg = projectNtjToAxesSmooth(  {nanmean(unitTimeTrialBZ,3)}, projVecs, gaussianKernel ); % project event(reward)-aligned trial-averaged trajectories
nTj.projToTagAxesAvgStmRch = projectNtjToAxesSmooth( {avgStmRchunitTimeTrialBZ}, projVecs, gaussianKernel ); % project stim reach trajectories
nTj.projToTagAxesAvgStmRwd = projectNtjToAxesSmooth( {avgStmRwdUnitTimeTrialBZ}, projVecs, gaussianKernel ); % project stim reward trajectories

%% visualize nTjs
% event markers and colorMaps
evtMarkers = [1, arrayfun(@(x) find(x==binEdges), evtMarkersRelative)]; % find the time bins corresponding to the evt markers using the relative evtMarkers (e.g. [-1050 0]) given as the input
nTjCmap = TNC_CreateRBColormapJP( p.Results.numbFolds, p.Results.trjCmap); 
nTjCmap = flip(nTjCmap, 1); % just to draw the faster/bigger reach trajectories with green  
evtCmap = TNC_CreateRBColormapJP( length(evtMarkers), 'cool');

tagEvtMarkers = arrayfun(@(x) find(x==tagRelativeTimeBins), tagWin); % find the time bins corresponding to the event markers for tagging activity trajectories
tagEvtCmap = repmat([255,165,0]./255, length(tagEvtMarkers),1); % just orange 

%% plot tagging activity + Trial-averaged reward-aligned + Trial-averaged reward-aligned laser On
figure; 
hold on; 
% plot the tagging activity trajectories projected onto the tagging axes
plot2DneuralTrajAndEventMarkers( {avgTagTrjMat}, [1,0,1], tagEvtMarkers, tagEvtCmap, p.Results.lineWidth, p.Results.markerSize ); 

% plot the event(reward)-aligned trajectories 'with light stimulation' projected onto the tagging axes
plot2DneuralTrajAndEventMarkers( nTj.projToTagAxesAvgStmRwd, [0 0 0], evtMarkers, evtCmap, p.Results.lineWidth, p.Results.markerSize ); 

% plot the event(reward)-aligned total trial-averaged trajectories
plot2DneuralTrajAndEventMarkers( nTj.projToTagAxesEvtAvg, [57 255 20]./255, evtMarkers, evtCmap, p.Results.lineWidth, p.Results.markerSize ); 

% plot the stmR trajectories projected onto the tagging axes
%stmRchbinEdges = StmRch.params.binEdges1ms(1):p.Results.binSize:StmRch.params.binEdges1ms(end)-1;
%stmRevtMarkers = [1, arrayfun(@(x) find(x==stmRchbinEdges), [-100 1450])]; % find the time bins corresponding to the evt markers using the relative evtMarkers (e.g. [-1050 0]) given as the input
% plot2DneuralTrajAndEventMarkers( nTj.projToTagAxesAvgStmRch, [0 0 0], stmRevtMarkers, evtCmap, p.Results.lineWidth, p.Results.markerSize ); 
hold off; 
print( fullfile(filePath,'Figure','tagAxesProjRwdAlignTrjs'), '-dpdf')

%% plot tagging activity + reward-aligned + velocity sorted/folded trjs
figure; 
hold on; 
% plot the tagging activity trajectories projected onto the tagging axes
plot2DneuralTrajAndEventMarkers( {avgTagTrjMat}, [1,0,1], tagEvtMarkers, tagEvtCmap, p.Results.lineWidth, p.Results.markerSize ); 

% plot the reward-aligned velocity sorted/folded trajectories
plot2DneuralTrajAndEventMarkers( nTj.projToTagAxesSortByMaxVel, nTjCmap, evtMarkers, evtCmap, p.Results.lineWidth, p.Results.markerSize ); 
title('TagDim_Tag+VelicityFolds','Interpreter', 'none')
hold off; 
print( fullfile(filePath,'Figure','tagAxesProjRwdAlignVelFoldTrjs'), '-dpdf')

%% plot tagging activity + reward-aligned + position sorted/folded trjs
figure; 
hold on; 
% plot the tagging activity trajectories projected onto the tagging axes
plot2DneuralTrajAndEventMarkers( {avgTagTrjMat}, [1,0,1], tagEvtMarkers, tagEvtCmap, p.Results.lineWidth, p.Results.markerSize ); 

% plot the reward-aligned velocity sorted/folded trajectories
plot2DneuralTrajAndEventMarkers( nTj.projToTagAxesSortByMaxPos, nTjCmap, evtMarkers, evtCmap, p.Results.lineWidth, p.Results.markerSize ); 
title('TagDim_Tag+PositionFolds','Interpreter', 'none')
hold off; 
print( fullfile(filePath,'Figure','tagAxesProjRwdAlignPosFoldTrjs'), '-dpdf')

%% plot tagging activity + reward-aligned + lickCount sorted/folded trjs
figure; 
hold on; 
% plot the tagging activity trajectories projected onto the tagging axes
plot2DneuralTrajAndEventMarkers( {avgTagTrjMat}, [1,0,1], tagEvtMarkers, tagEvtCmap, p.Results.lineWidth, p.Results.markerSize ); 

% plot the reward-aligned velocity sorted/folded trajectories
plot2DneuralTrajAndEventMarkers( nTj.projToTagAxesSortByLickCnt, nTjCmap, evtMarkers, evtCmap, p.Results.lineWidth, p.Results.markerSize ); 
title('TagDim_Tag+LickCountFolds','Interpreter', 'none')
hold off; 
print( fullfile(filePath,'Figure','tagAxesProjRwdAlignLickCntFoldTrjs'), '-dpdf')

%% save data
save(fullfile(filePath,strcat(saveNameTag,'_projectTo2dTaggingAxes')), 'nTj', 'bTj', 'avg*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_projectTrialAveragedNtjTo2dTaggingAxis( filePath, tagPCAFileName, binSpkCntFileName, evtName, saveNameTag, vargs )
        
        default_binSize = 50; % 50 ms bins
        default_kernelSD = 100; % 100 ms sd for the gaussian kernel used for smoothing
        default_numbFolds = 3;   % # of folds to divide nTjs into
        default_lineWidth  = 1;  % default lineWidth to be used to plot the neural trajectories
        default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
        default_trjCmap = 'summer'; % the default trajectory colormap
        default_PCs = 1:2; % PC directions (axes) to project nTrjs onto
        
        p = inputParser; % create parser object
        
        addRequired(p,'filePath'); % file directory
        addRequired(p,'tagPCAFileName'); % fileName for tagging acticity pca result
        addRequired(p,'binSpkCntFileName'); % fileName for binSpkCnt
        addRequired(p,'saveNameTag'); % saveName used to save the outcomes 
        addRequired(p,'evtName'); % event name the neural trajectories are aligned to
        
        addParameter(p,'binSize', default_binSize) 
        addParameter(p,'kernelSD', default_kernelSD)
        addParameter(p,'numbFolds', default_numbFolds) 
        addParameter(p,'lineWidth', default_lineWidth)   
        addParameter(p,'markerSize', default_markerSize) 
        addParameter(p,'trjCmap', default_trjCmap) 
        addParameter(p,'PCs', default_PCs) 
        
        parse(p, filePath, tagPCAFileName, binSpkCntFileName, evtName, saveNameTag, vargs{:})
    end

end

