%This script conducts a targeted dimensionality reduction analysis to
% define the direction within a population state space, along which the
% neural population activity covaries as a function of movement kinematic
% variables.
filePath = '/Volumes/dudmanlab/junchol/data_backup/ITPTphys/ITPhys/IT01_121317/Matfiles'; 
cd(filePath)

%unitIdx = load(fullfile(filePath,'meanFRxCorrUnitIdxCtx.mat')); % load unitIdx  
%unitIdx = unitIdx.('meanFRxCorrUnitIdx'); % units to be used 
fileNameUnitTimeTrial = 'binSpkCountCTX_KS_IT01_121317.mat'; 
fileNameNeuralTrj     = 'IT01Ctx_121317_pca_reach_baseSub_KS_5D_50msBin.mat';
fileNameTempPCA       = 'binSpkCountStrCtx_KS_IT01_121317pcaPSTHreach200ms.mat';
fileNameNeuralTagTrj  = 'IT01Ctx_121317_pca_tag_baseSub_KS_5D_50msBin.mat'; 
saveName = 'IT01Ctx_121317_targetedDimReductionRegression_reachAligned_stmLaser'; 

trialFolds = 3; % # of folds to divide the kinematic variables with

%% determine ETag and ITag dimensions
% this is just to see which tagging dimension captures the excitatory/inhibitory tagging effects
[nTrjCtxTag] = visualizeNeuralTraj(filePath, 'IT01Ctx_121317_pca_tag_baseSub_KS_5D_50msBin.mat', 'IT01Ctx_121317_pca_tag_baseSub_KS_5D_50msBin', [0 1000] ...
,'whatToPlot','Folds','numbTrials', 10,'dimPlotBy', 2,'numbFolds',1,'trjCmap','summer'); 
nTrjCtxTagSortedFold = nTrjCtxTag.nTrjCell; 
hold on; 
plotNtrjEachDimByTime(filePath, nTrjCtxTagSortedFold, 1, nTrjCtxTag.relativeTimeBins, 'summer') % plot PCscores of each dimension over time (dim #1)
plotNtrjEachDimByTime(filePath, nTrjCtxTagSortedFold, 2, nTrjCtxTag.relativeTimeBins, 'summer') % plot PCscores of each dimension over time (dim #2)
hold off; 

eTagD = 2;
iTagD = 1;

%% load files 
% get behavioral data 'BehVariables.mat'
load(fullfile(filePath,'BehVariables.mat'),'ts','reach0','lick')
ts.reachNoStimIdx = ~logical(cellfun(@(x) sum(abs(x-ts.stmLaser)<=2000), num2cell(ts.reachStart))); % index for reward trials stim off

% get 
S = load(fullfile(filePath,fileNameTempPCA),'S'); 
S = S.('S'); 
if isfield(S,'nonNaNtrialId')
    nonNaNtrID = S.nonNaNtrialId;  
end
clearvars S 

% get the unitTimeTrial mat for reward aligned neural activity
S = load(fullfile(filePath,fileNameUnitTimeTrial),'reach'); 
S = S.('reach'); 

unitTimeTrialRch = S.unitTimeTrial; %clearvars S, REPLACE unitTimeTrialRwd to unitTimeTrialRch

% get the unitTimeTrial mat for stmReach aligned neural activity
S = load(fullfile(filePath,fileNameUnitTimeTrial),'stmLaser'); 
S = S.('stmLaser'); 
unitTimeTrialStmLaser = S.unitTimeTrial; clearvars S

% get pca- or gpfa-based neural trajectories 
pcaResultNTj = load(fileNameNeuralTrj,'pcaResult'); 
pcaResultNTj = pcaResultNTj.('pcaResult'); 

pcaResultTag = load(fileNameNeuralTagTrj,'pcaResult'); 
pcaResultTag = pcaResultTag.('pcaResult'); 
unitIdx = pcaResultTag.unitIdx; % use the pcaResultTag.unitIdx to match the dimension of the tagging activity PCs (to project onto them)

% sanity check for trialId
if size(unitTimeTrialRch,3)~=length([pcaResultNTj.kern.seqTrain.trialId]) % trial IDs for the neural trjectories
   if length([pcaResultNTj.kern.seqTrain.trialId]) == length(nonNaNtrID)
       trialId = [pcaResultNTj.kern.seqTrain.trialId]; 
       unitTimeTrialRch = unitTimeTrialRch(:,:,nonNaNtrID); 
   else
       error('The # of trials do not match!')
   end
elseif size(unitTimeTrialRch,3)==length([pcaResultNTj.kern.seqTrain.trialId])
    trialId = [pcaResultNTj.kern.seqTrain.trialId]; 
end

% load stimE
load(fullfile(filePath, 'binSpkCountCTX_KS_IT01_121317_stimE.mat'))
if ~(sum(stimE.tagE(unitIdx)==1)>=10 && sum(stimE.tagE(unitIdx)==-1)>=10)
    warning('There is not enough number of tagged units to comprise tagE and tagI dimensions!')
end

%% get the event windows
if contains(fileNameNeuralTrj,'reach') && ~contains(fileNameNeuralTrj,'reward')
    behTS = ts.reachStart(trialId); % trialId (*seqTrain.trialId) contains non-NaN trials only
    eventMarkersRelative = [-50 1450]; % reachStart, reward delivery (roughly) 
    reachWin = [-200 1000]; % reach window relative to ReachStart
    lickWin = [500 3000]; % reward window relatie to ReachStart
elseif contains(fileNameNeuralTrj,'reward') && ~contains(fileNameNeuralTrj,'reach')
    behTS = ts.reward(trialId); % trialId (*seqTrain.trialId) contains non-NaN trials only
    eventMarkersRelative = [-1600 0]; 
    reachWin = [-1500 0]; % reach window relative to Reward
    lickWin = [0 2500]; % reward window relative to Reward
else
    error('Make sure if correct neural data has been input!')
end

%% preprocess unitTimeTrial 
% bin the spike count mat with 50-ms sliding window
unitTimeTrialB = bin1msSpkCountMat( unitTimeTrialRch, 50, 50, 'align', 'center' ); 

% get the mean and std across all 50-ms time bins and trials, then z-score normalize
rsUnitTimeTrial = reshape(unitTimeTrialB, [size(unitTimeTrialB,1), size(unitTimeTrialB,2)*size(unitTimeTrialB,3)]); 
meanPerUnit = nanmean(rsUnitTimeTrial,2); 
stdPerUnit = nanstd(rsUnitTimeTrial,0,2); 
rsUnitTimeTrialBZ=(rsUnitTimeTrial-repmat(meanPerUnit,[1 size(rsUnitTimeTrial,2)]))./repmat(stdPerUnit,[1 size(rsUnitTimeTrial,2)]); 
unitTimeTrialBZ = reshape(rsUnitTimeTrialBZ, [size(unitTimeTrialB,1), size(unitTimeTrialB,2), size(unitTimeTrialB,3)]); 
%unitTimeTrialBZ = unitTimeTrialBZ(:,:,ts.reachNoStimIdx); 

%% preprocess (non-targeted) pca- or gpfa-based neural trajectories 
nTj.trjMat  = reshape([pcaResultNTj.kern.seqTrain.xpost], pcaResultNTj.p.Results.pcaDim, [], length(pcaResultNTj.kern.seqTrain)); % neural trajectories in dim x timeBins x trials (e.g. 5x100x114)
nTj.trialId = [pcaResultNTj.kern.seqTrain.trialId]; % trial IDs for the neural trjectories
nTj.relativeTimeBins = pcaResultNTj.p.Results.timeRange(1):pcaResultNTj.p.Results.pcaBinSize:pcaResultNTj.p.Results.timeRange(2)-pcaResultNTj.p.Results.pcaBinSize; % get timeBins relative to time 0 (e.g. reward) of the current PSTH
nTj.reachTimeBins = nTj.relativeTimeBins>=reachWin(1) & nTj.relativeTimeBins<reachWin(2); % the time bins of interest for the behavioral kinematics (e.g. when reaches are most likely to occur) to narrow down the range of behavioral measures
nTj.lickTimeBins  = nTj.relativeTimeBins>=lickWin(1) & nTj.relativeTimeBins<lickWin(2); % the time bins of interest for the behavioral kinematics (e.g. when reaches are most likely to occur) to narrow down the range of behavioral measures
nTj.neuralTrajFile = fullfile(filePath,fileNameNeuralTrj); % to keep track of the neural trajectory files
nTj.pcaRez = pcaResultNTj.kern; % pca info pc loadings, eigVals, expVar 

%% preprocess the behavioral data
lickBin = zeros(1,length(lick)); % licks
lickBin(ts.lick) = 1; % licks

valTrialIdx = zeros(length(behTS),1);
for t = 1:length(behTS) % increment trials, take the trial-by-trial position/velocity data and bin them to match the neural trjectories (e.g. 50 ms)
    if max(behTS(t)+pcaResultNTj.p.Results.timeRange)<=length(reach0) && min(behTS(t)+pcaResultNTj.p.Results.timeRange)>0
        valTrialIdx(t)=1;
        timeWin = behTS(t)+pcaResultNTj.p.Results.timeRange; % the time window, e.g. -3 to 2 sec relative to the behavioral timestamp
        bTj(t).reachPos = smooth(decimate(reach0(timeWin(1):timeWin(2)-1),pcaResultNTj.p.Results.pcaBinSize),3)'; % get decimated behavioral trjectories on the same timescale of the neural trjectories
        bTj(t).reachVel = smooth(diff([bTj(t).reachPos(1) bTj(t).reachPos]),3)'.*1000; % get reach velocities (from reach0)
        bTj(t).maxReachPos = max(bTj(t).reachPos(nTj.reachTimeBins)); % max reach position within the time window of interest
        [~,bTj(t).maxReachPosI] = max(bTj(t).reachPos(nTj.reachTimeBins)); % max reach position bin Idx 
        bTj(t).maxReachVel = max(bTj(t).reachVel(nTj.reachTimeBins)); % max reach velocity within the time window of interest
        bTj(t).trialId = t;
        bTj(t).lick = bin1msSpkCountMat(lickBin(timeWin(1):timeWin(2)-1),pcaResultNTj.p.Results.pcaBinSize,pcaResultNTj.p.Results.pcaBinSize); % get binned lick counts using the bin1msSpkCountMat
        bTj(t).lickCount =  sum(bTj(t).lick(nTj.lickTimeBins)); % lick Counts within the lickTimeBins
    else
    end
end
clearvars t
maxReachPosIcrt = find(-3000:50:2000 == reachWin(1))+[bTj.maxReachPosI]-1; % binID for maxReachPosition 

%% get regression coefficients per unit per time bin (Bit, i; unitIdx, t; timeBinIdx)
% get regressors (position, velocity, lick count)
maxP = [[bTj.maxReachPos];[bTj.trialId]]'; % trial-by-trial reach positions
maxV = [[bTj.maxReachVel];[bTj.trialId]]'; % trial-by-trial max reach velocities
lickCnt = [[bTj.lickCount];[bTj.trialId]]'; % trial-by-trial lick counts

maxP(:,3) = discretize(maxP(:,1),linspace(min(maxP(:,1)),max(maxP(:,1)),trialFolds+1)); % discretize trial-by-trial positions
maxV(:,3) = discretize(maxV(:,1),linspace(min(maxV(:,1)),max(maxV(:,1)),trialFolds+1)); % discretize trial-by-trial velocities
lickCnt(:,3) = discretize(lickCnt(:,1),linspace(min(lickCnt(:,1)),max(lickCnt(:,1)),trialFolds+1)); % discretize trial-by-trial lick counts

FiT = maxP(:,3); % discretized regressors
%FiT = [maxP(:,3),maxV(:,3),lickCnt(:,3)]; % discretized regressors
FiT(:,end+1) = ones(size(FiT,1),1); % add ones in the last row
Fi = FiT'; % just transpose the regressors

invFftf = (Fi*Fi')\Fi; % inv(Fi*Fi')*Fi; 
permUnitTimeTrialBZ = permute(unitTimeTrialBZ(unitIdx,:,:),[3,1,2]); 
rsPermUnitTimeTrialBZ = reshape(permUnitTimeTrialBZ, [size(permUnitTimeTrialBZ,1),size(permUnitTimeTrialBZ,2)*size(permUnitTimeTrialBZ,3)]);
Bit = reshape(invFftf*rsPermUnitTimeTrialBZ, [size(Fi,1), sum(unitIdx), size(unitTimeTrialBZ,2)]); 
% Bit(:,:,4)-invFftf*squeeze(permUnitTimeTrialB(:,:,4)) % to verify the regression coefficient calculation 

maxBi = max(abs(Bit),[],3); % dimension (#Regressors-by-#Units) argmax regressors across time bins

bRg.maxP = maxP; 
bRg.maxV = maxV; % put variables relevant to regression analysis into bTjR  
bRg.lickCnt = lickCnt; 
bRg.Fi = Fi; 
bRg.Bit = Bit; 
bRg.maxBi = maxBi; 

%% use Matlab's regress to do a significance test on the regression coefficients
% get the unitTimeTrial mat for reach aligned neural activity
S = load(fullfile(filePath,fileNameUnitTimeTrial),'reach'); 
S = S.('reach'); 
trAvgUnitTimeTrialBZreach = cell2mat(cellfun(@(a) binAvg1msSpkCountMat(a, 50, 50), S.SpkCountMatZ, 'Un', 0));
trAvgUnitTimeTrialBZreach = trAvgUnitTimeTrialBZreach(unitIdx,:); 

permUnitTimeTrialBZLmt = permUnitTimeTrialBZ(:,:,10:80); 
for u = 1:size(permUnitTimeTrialBZLmt, 2) % increment unit
    for t = 1:size(permUnitTimeTrialBZLmt, 3) % increment time
        [rcoeff,~,~,~,rStats] = regress(permUnitTimeTrialBZLmt(valTrialIdx==1,u,t), FiT); 
        reg.B(u,t)= rcoeff(1); 
        reg.rsq(u,t) = rStats(1); 
        reg.Fval(u,t) = rStats(2); 
        reg.Pval(u,t) = rStats(3); 
        % corr
        [cor.rho(u,t),cor.pval(u,t)] = corr(permUnitTimeTrialBZLmt(valTrialIdx==1,u,t),FiT(:,1)); 
    end
    [~,reg.maxTI(u,1)] = max(abs(reg.B(u,:)),[],2);
    reg.absMaxB(u,1) = reg.B(u,reg.maxTI(u,1)); 
    reg.absMaxP(u,1) = reg.Pval(u,reg.maxTI(u,1)); 
    reg.absMaxRsq(u,1) = reg.rsq(u,reg.maxTI(u,1)); 
    reg.zMove(u,1) = nanmean(trAvgUnitTimeTrialBZreach(u,30:60))-nanmean(trAvgUnitTimeTrialBZreach(u,1:40)); 
    [~,reg.zMoveAbsMaxI(u,1)] = max(abs(trAvgUnitTimeTrialBZreach(u,40:60))); 
    reg.zMoveAbsMax(u,1) = trAvgUnitTimeTrialBZreach(u,reg.zMoveAbsMaxI(u,1)); 
    % run regression with the summed spike count in the movement window
    [rcoeffB,~,~,~,rStatsB] = regress(sum(squeeze(unitTimeTrialB(u,20:60,valTrialIdx==1)),1)', FiT); 
    reg.moveBinB(u,1) = rcoeffB(1); 
    reg.moveBinP(u,1) = rStatsB(3); 
    
    %sigCor = cor.rho(u,cor.pval(u,:)<.05); 
    %[~, cor.absMaxRI(u,1)] = max(abs(sigCor),[],2); 

end
clearvars u 

reg.itI = stimE.tagE(unitIdx)==-1; 
reg.etI = stimE.tagE(unitIdx)==1; 
reg.utI = stimE.tagE(unitIdx)==0; 
reg.sigI = reg.absMaxP<0.05; 
reg.sigIbin = reg.moveBinB<.05; 

randX = -.1 + (.1+.1)*rand(1000,1); 
% plot tag effect
scatterColorSpace = linspace(-5,5,20); % to color code each scatter by the tagEHalfPeakT
scatterColor = nan(length(reg.moveBinB),1);
for u = 1:length(reg.moveBinB)
    if reg.zMove(u)<=min(scatterColorSpace)
      scatterColor(u,1) = scatterColorSpace(1);   
    else
      scatterColor(u,1) = scatterColorSpace(find(scatterColorSpace<=reg.zMove(u),1,'last')); % match the color from the scatterColorSpace
    end
end

figure; hold on; 
scatter(1+randX(reg.itI&reg.sigI),reg.absMaxB(reg.itI&reg.sigI),50,scatterColor(reg.itI&reg.sigI),'filled');
scatter(1+randX(reg.itI&~reg.sigI),reg.absMaxB(reg.itI&~reg.sigI),50,scatterColor(reg.itI&~reg.sigI));

scatter(2+randX(reg.utI&reg.sigI),reg.absMaxB(reg.utI&reg.sigI),50,scatterColor(reg.utI&reg.sigI),'filled');
scatter(2+randX(reg.utI&~reg.sigI),reg.absMaxB(reg.utI&~reg.sigI),50,scatterColor(reg.utI&~reg.sigI));
colormap jet
caxis([-5 5])

% figure; hold on; 
% scatter(1+randX(reg.ptI&reg.sigIbin),reg.moveBinB(reg.ptI&reg.sigIbin),50,scatterColor(reg.ptI&reg.sigIbin),'filled');
% scatter(1+randX(reg.ptI&~reg.sigIbin),reg.moveBinB(reg.ptI&~reg.sigIbin),50,scatterColor(reg.ptI&~reg.sigIbin));
% 
% scatter(2+randX(reg.itI&reg.sigIbin),reg.moveBinB(reg.itI&reg.sigIbin),50,scatterColor(reg.itI&reg.sigIbin),'filled');
% scatter(2+randX(reg.itI&~reg.sigIbin),reg.moveBinB(reg.itI&~reg.sigIbin),50,scatterColor(reg.itI&~reg.sigIbin));
% colormap jet
% hold off; 

%% population average responses
% get the population average responses of all combinations of behavioral variables
popAvgC = cell(trialFolds,trialFolds,trialFolds); 
pcaPopAvg = []; % take population 

% gaussian kernel to be convolved with the psths
gaussianSigma    = 1;  % gaussian std (50ms, as the data are already binned with 50-ms window)
gaussianKernel  = TNC_CreateGaussian(gaussianSigma*15,gaussianSigma,gaussianSigma*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)
gaussianKernel2 = TNC_CreateGaussian(gaussianSigma*15,gaussianSigma*2,gaussianSigma*30,1); 

% get the unit-by-timeBin population activity mat averaged across trials
%within each condition (every combination of discretized position-velocity-lickCnt values)
for p = 1:trialFolds 
    tmpP = p; 
    for v = 1:trialFolds
        tmpV = v;
        for l = 1:trialFolds
            tmpL = l;
            popAvgC{p,v,l} = unitTimeTrialBZ(unitIdx,:,maxP(:,3)==tmpP & maxV(:,3)==tmpV & lickCnt(:,3)==tmpL);
            if ~isempty(popAvgC{p,v,l})
                tmpPopAvg = nanmean(popAvgC{p,v,l},3); % average across trials within each condition
                smTmpPopAvg = cell2mat(arrayfun(@(ROWIDX) conv(tmpPopAvg(ROWIDX,:),gaussianKernel,'same'), (1:size(tmpPopAvg,1)).', 'UniformOutput', false)); % smoothing by convolving with the gaussian kernel
                pcaPopAvg = [pcaPopAvg, smTmpPopAvg]; % append smTmpPopAvg to construct the pcaPopAvg (unit x (timeBin x conditions))
            end
        end 
    end
end
clearvars p v l

%% targeted dimensionality reduction  
% 1) denoising the regression vectors (coefficients) by projecting onto the subspace spanned by the first 10 PCs. 
[popAvgPCdir] = pca(pcaPopAvg'); % unit-by-unit pc dir
% construct the denoising mat using the first 10 PCs
denoisePCdir = zeros(size(popAvgPCdir,1), size(popAvgPCdir,2),10); 
for pc = 1:10
    denoisePCdir(:,:,pc) = popAvgPCdir(:,pc)*popAvgPCdir(:,pc)'; 
end
clearvars pc
D = sum(denoisePCdir,3); % the denoising matrix

Bvt = permute(Bit,[2 1 3]); % just rearranging the Bit to Bvt which is in #Units-by-#Regressors-by-#TimeBins
for t = 1:size(Bvt,3) % increment timeBins
    BvtD(:,:,t) = D*Bvt(:,:,t); % denoised regression coefficeints 
end

% 2) get the time-independent denoised regression vectors
maxBvtD = max(BvtD,[],3); % the time-independent, denoised regression vectors

% 3) orthogonalize (to ensure that each axis explains distinct portions of the variance) - QR decomposition
[Q,R] = qr(maxBvtD); 
Bvortho = Q(:,1:size(maxBvtD,2)); 

%% project to tagging activity PCs
projectMat = [pcaResultTag.kern.pcDirs(:,[eTagD,iTagD]), -Bvortho(:,1)]; % construct the axes to which neural population activity will be projected (e.g. velocity axis, first two PCs of the tagging activity)
for f = 1:trialFolds
    tmpNtj = projectMat'*nanmean(unitTimeTrialBZ(unitIdx,:,maxP(:,3)==f),3); % project the trial-averaged population activity of each fold to the projection matrix
    nTjCell{1,f} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtj(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtj,1))', 'UniformOutput', false)); % smooth
end
clearvars f

%% visualize trajectories
eventMarkers = [1, arrayfun(@(x) find(x==nTj.relativeTimeBins), eventMarkersRelative)]; % time points to mark on the neural trajectories
nTjCmap = flip(TNC_CreateRBColormapJP( length(nTjCell), 'summer'),1); % flip to keep up with the convention that darker green corresponds to greater kinematic values
nTjCmap2 = [[39, 170, 225]./255; [79, 138, 205]./255; [21, 67, 130]./255]; 
nTjCmap3 = [[65, 64, 66]./255; [158, 31, 99]./255]; 

evtCmap = TNC_CreateRBColormapJP( length(eventMarkers), 'cool'); 
plot3DneuralTrajAndEventMarkers( nTjCell, nTjCmap, eventMarkers, evtCmap, 2, 10 )

for tf = 1:trialFolds
    projP{tf} = -conv(Bvortho(:,1)'*nanmean(unitTimeTrialBZ(unitIdx,:,maxP(:,3)==tf),3),gaussianKernel,'same');  
    projETag{tf} = conv(projectMat(:,1)'*nanmean(unitTimeTrialBZ(unitIdx,:,maxP(:,3)==tf),3),gaussianKernel,'same');
    projITag{tf} = conv(projectMat(:,2)'*nanmean(unitTimeTrialBZ(unitIdx,:,maxP(:,3)==tf),3),gaussianKernel,'same');
    projTaskPC{1,tf} = conv(pcaResultNTj.kern.pcDirs(:,1)'*nanmean(unitTimeTrialBZ(pcaResultNTj.unitIdx,:,maxP(:,3)==tf),3),gaussianKernel,'same'); 
    projTaskPC{2,tf} = conv(pcaResultNTj.kern.pcDirs(:,2)'*nanmean(unitTimeTrialBZ(pcaResultNTj.unitIdx,:,maxP(:,3)==tf),3),gaussianKernel,'same'); 
    projTaskPC{3,tf} = conv(pcaResultNTj.kern.pcDirs(:,3)'*nanmean(unitTimeTrialBZ(pcaResultNTj.unitIdx,:,maxP(:,3)==tf),3),gaussianKernel,'same');    
    c{tf} = sprintf('trFold#%d',tf);  
end
clearvars tf

xAxis = nTj.relativeTimeBins; 
clear g
g(1,1)=gramm('x',xAxis,'y',projP,'color',c);
g(1,1).geom_line(); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,1).set_names('x','Time (ms)','y','Projection to Dist. Dim a.u.');
g(1,1).set_title('Distance');  

g(1,2)=gramm('x',xAxis,'y',projETag,'color',c);
g(1,2).geom_line(); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,2).set_names('x','Time (ms)','y','Projection to ETag Dim a.u.');
g(1,2).set_title('ETag');  
 
g(1,3)=gramm('x',xAxis,'y',projITag,'color',c);
g(1,3).geom_line(); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,3).set_names('x','Time (ms)','y','Projection to ITag Dim a.u.');
g(1,3).set_title('ITag');   

figHandle = figure('Position',[100 100 800 350]);
g.set_title('Projection to targeted kinematic dimensions'); 
g.draw();

print( fullfile(filePath,'Figure',strcat('projectionToTargetedKinematicDimensions')), '-dpdf', '-bestfit')

%% save info neural activity projected to each dimensions
dotProdExtTagDistanceDim = abs(pcaResultTag.kern.pcDirs(:,eTagD)'*Bvortho(:,1)); 
dotProdInhTagDistanceDim = abs(pcaResultTag.kern.pcDirs(:,iTagD)'*Bvortho(:,1)); 

maxRegCoeff.NegTag = maxBi(1,stimE.tagE(unitIdx)==-1); 
maxRegCoeff.PosTag = maxBi(1,stimE.tagE(unitIdx)==1); 
maxRegCoeff.NotTag = maxBi(1,stimE.tagE(unitIdx)==0); 

save(fullfile(filePath,saveName))

%% project to tagging activity PCs and movement dimension
for t = 1:size(unitTimeTrialBZ,3) % trial-by-trial
    tmpNtj = projectMat'*unitTimeTrialBZ(unitIdx,:,t); % project the trial-averaged population activity of each fold to the projection matrix
    nTjCellAllTrials{1,1,t} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtj(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtj,1))', 'UniformOutput', false)); % smooth
end
clearvars t

valStimTrI = ts.reachNoStimIdx(trialId); 

prjKN_nPtb = cell2mat(squeeze(cellfun(@(a) a(3,:), nTjCellAllTrials(valStimTrI), 'un', 0))); 
prjKN_ptb = cell2mat(squeeze(cellfun(@(a) a(3,:), nTjCellAllTrials(~valStimTrI), 'un', 0))); 

nTjC.nPtb = nanmean(prjKN_nPtb);
nTjC.ptb = nanmean(prjKN_ptb);

nTjCellAllTrialsExtMovDim = cellfun(@(x) x([1,3],:), nTjCellAllTrials, 'UniformOutput', false); % project onto ExTag x Distance dims 
nTjCellAllTrialsInhMovDim = cellfun(@(x) x([2,3],:), nTjCellAllTrials, 'UniformOutput', false); % project onto InhTag x Distance dims

nTjCellAllTrialsExtMovDimMeanNperturb = nanmean(cell2mat(nTjCellAllTrialsExtMovDim(1,1,valStimTrI)),3); 
nTjCellAllTrialsExtMovDimMeanPerturb = nanmean(cell2mat(nTjCellAllTrialsExtMovDim(1,1,~valStimTrI)),3); 

nTjCellAllTrialsInhMovDimMeanNperturb = nanmean(cell2mat([nTjCellAllTrialsInhMovDim(1,1,valStimTrI)]),3); 
nTjCellAllTrialsInhMovDimMeanPerturb = nanmean(cell2mat([nTjCellAllTrialsInhMovDim(1,1,~valStimTrI)]),3); 

evtMarkersRelative = [-1000 0]-50; % time points to mark on the neural trajectories
evtMarkers = [1, arrayfun(@(x) find(x==[-3000:50:2000]), evtMarkersRelative)]; % find the time bins corresponding to the evt markers using the relative evtMarkers (e.g. [-1050 0]) given as the input
%nTjCmap = []./255; 
evtCmap = TNC_CreateRBColormapJP( length(evtMarkers), 'cool');

%% project unitTimeTrialBZ_stmLaser to tagging activity and KN dims
% get the mean and std across all 50-ms time bins and trials, then z-score normalize
unitTimeTrialB_stmLaser = bin1msSpkCountMat( unitTimeTrialStmLaser, 50, 50, 'align', 'center' ); 
rsUnitTimeTrial_stmLaser = reshape(unitTimeTrialB_stmLaser, [size(unitTimeTrialB_stmLaser,1), size(unitTimeTrialB_stmLaser,2)*size(unitTimeTrialB_stmLaser,3)]); 
meanPerUnit = nanmean(rsUnitTimeTrial_stmLaser,2); 
stdPerUnit = nanstd(rsUnitTimeTrial_stmLaser,0,2); 
rsUnitTimeTrial_stmLaserB_stmLaserZ=(rsUnitTimeTrial_stmLaser-repmat(meanPerUnit,[1 size(rsUnitTimeTrial_stmLaser,2)]))./repmat(stdPerUnit,[1 size(rsUnitTimeTrial_stmLaser,2)]); 
unitTimeTrialBZ_stmLaser = reshape(rsUnitTimeTrial_stmLaserB_stmLaserZ, [size(unitTimeTrialB_stmLaser,1), size(unitTimeTrialB_stmLaser,2), size(unitTimeTrialB_stmLaser,3)]); 

for t = 1:size(unitTimeTrialBZ_stmLaser,3) % trial-by-trial
    tmpNtj_stmLaser = projectMat'*unitTimeTrialBZ_stmLaser(unitIdx,:,t); % project the trial-averaged population activity of each fold to the projection matrix
    nTjCellAllTrials_stmLaser{1,1,t} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtj_stmLaser(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtj_stmLaser,1))', 'UniformOutput', false)); % smooth
end
clearvars t

prjKN_stmLaser = cell2mat(squeeze(cellfun(@(a) a(3,:), nTjCellAllTrials_stmLaser, 'un', 0))); 
nTjC.stmLaser = nanmean(prjKN_stmLaser);

save(fullfile(filePath,saveName), 'nTjC', 'prjKN*', '-append')

%% plot neural trajectories projected onto each tag dim + distance dim
nTjNptbPtbC = {nTjC.nPtb(:,1:60),nTjC.ptb(:,1:60)}; 
figure; 
plot3DneuralTrajAndEventMarkers( nTjNptbPtbC, nTjCmap3, evtMarkers, evtCmap, 1.5, 8 ); 
% threshold crossing points on each axis
hold on; 
plot3(ones(1,40)*nTjNptbPtbC{1,1}(1,40), ones(1,40)*nTjNptbPtbC{1,1}(2,40),-2:(nTjNptbPtbC{1,1}(3,40)+2)/40:nTjNptbPtbC{1,1}(3,40)-(nTjNptbPtbC{1,1}(3,40)+2)/40, 'Color',[65, 64, 66]./255,'LineStyle',':','LineWidth',1); 
plot3(ones(1,40)*nTjNptbPtbC{1,2}(1,40), ones(1,40)*nTjNptbPtbC{1,2}(2,40),-2:(nTjNptbPtbC{1,2}(3,40)+2)/40:nTjNptbPtbC{1,2}(3,40)-(nTjNptbPtbC{1,2}(3,40)+2)/40, 'Color',[158, 31, 99]./255,'LineStyle',':','LineWidth',1); 

plot3(-2:(nTjNptbPtbC{1,1}(1,40)+2)/40:nTjNptbPtbC{1,1}(1,40)-(nTjNptbPtbC{1,1}(1,40)+2)/40, ones(1,40)*nTjNptbPtbC{1,1}(2,40),ones(1,40)*-2, 'Color',[65, 64, 66]./255,'LineStyle',':','LineWidth',1); 
plot3(-2:(nTjNptbPtbC{1,2}(1,40)+2)/40:nTjNptbPtbC{1,2}(1,40)-(nTjNptbPtbC{1,2}(1,40)+2)/40, ones(1,40)*nTjNptbPtbC{1,2}(2,40),ones(1,40)*-2, 'Color',[158, 31, 99]./255,'LineStyle',':','LineWidth',1); 

plot3(ones(1,40)*nTjNptbPtbC{1,1}(1,40), -3:(nTjNptbPtbC{1,1}(2,40)+3)/40:nTjNptbPtbC{1,1}(2,40)-(nTjNptbPtbC{1,1}(2,40)+3)/40,ones(1,40)*-2, 'Color',[65, 64, 66]./255,'LineStyle',':','LineWidth',1); 
plot3(ones(1,40)*nTjNptbPtbC{1,2}(1,40), -3:(nTjNptbPtbC{1,2}(2,40)+3)/40:nTjNptbPtbC{1,2}(2,40)-(nTjNptbPtbC{1,2}(2,40)+3)/40,ones(1,40)*-2, 'Color',[158, 31, 99]./255,'LineStyle',':','LineWidth',1); 
print( fullfile(filePath,'Figure',strcat('projectionToETag-ITag-DistDimensions_perturbUnperturb')), '-dpdf', '-bestfit', '-painters')

figure; 
hold on; % Ext-Dim x Distance-Dim 
% plot the perturbed vs. non-perturbed trials (trial-averaged)
plot2DneuralTrajAndEventMarkers( {nTjCellAllTrialsExtMovDimMeanNperturb}, [0,1,1], evtMarkers, evtCmap, 1, 8 ); 
plot2DneuralTrajAndEventMarkers( {nTjCellAllTrialsExtMovDimMeanPerturb}, [1,0,1], evtMarkers, evtCmap, 1, 8 ); 
xlim([-2 4])
ylim([-2 4])
hold off; 
print( fullfile(filePath,'Figure',strcat('projectionToETag-DistDimensions_Perturb')), '-dpdf', '-bestfit')

figure; % Inh-Dim x Distance-Dim
hold on; 
% plot the perturbed vs. non-perturbed trials (trial-averaged)
plot2DneuralTrajAndEventMarkers( {nTjCellAllTrialsInhMovDimMeanNperturb}, [0,1,1], evtMarkers, evtCmap, 1, 5 ); 
plot2DneuralTrajAndEventMarkers( {nTjCellAllTrialsInhMovDimMeanPerturb}, [1,0,1], evtMarkers, evtCmap, 1, 5 ); 
xlim([-4 2]);
ylim([-2 4]);
hold off; 
print( fullfile(filePath,'Figure',strcat('projectionToITag-DistDimensions_Perturb')), '-dpdf', '-bestfit')

figure; % ETag-Dim x Distance-Dim Distance tertiles
hold on; 
% plot the max distance tertiles (trial-averaged)
plot2DneuralTrajAndEventMarkers( {nTjCell{1}([1,3],:)}, [39, 170, 225]./255, evtMarkers, evtCmap, 1, 5 ); 
plot2DneuralTrajAndEventMarkers( {nTjCell{2}([1,3],:)}, [79, 138, 205]./255, evtMarkers, evtCmap, 1, 5 ); 
plot2DneuralTrajAndEventMarkers( {nTjCell{3}([1,3],:)}, [21, 67, 130]./255, evtMarkers, evtCmap, 1, 5 ); 
xlim([-2 4]);
ylim([-2 4]);
hold off; 
print( fullfile(filePath,'Figure',strcat('projectionToETag-DistDimensions')), '-dpdf', '-bestfit')


nTjCellCut = cellfun(@(c) c(:,1:60), nTjCell,'Un',0); % upto the reward delivery point
figure; % ETag-ITag-Dim x Distance-Dim Distance tertiles
plot3DneuralTrajAndEventMarkers( nTjCellCut, nTjCmap2, evtMarkers, evtCmap, 1.5, 8 ); 
maxVal = cellfun(@(a) max(a,[],2), nTjCell, 'Un', 0);
maxValMat = max(cell2mat(maxVal),[],2); 
minVal = cellfun(@(a) min(a,[],2), nTjCell, 'Un', 0);
minValMat = min(cell2mat(minVal),[],2);
%pbaspect(maxValMat-minValMat); 
pbaspect([6 2 6]); 
valAtT = cell2mat(cellfun(@(a) a(:,40), nTjCell, 'Un', 0)); 
hold on; 
% threshold crossing points on each axis
plot3(ones(1,40)*nTjCellCut{1,1}(1,40), ones(1,40)*nTjCellCut{1,1}(2,40),-2:(nTjCellCut{1,1}(3,40)+2)/40:nTjCellCut{1,1}(3,40)-(nTjCellCut{1,1}(3,40)+2)/40, 'Color',[39, 170, 225]./255,'LineStyle',':','LineWidth',1); 
plot3(ones(1,40)*nTjCellCut{1,2}(1,40), ones(1,40)*nTjCellCut{1,2}(2,40),-2:(nTjCellCut{1,2}(3,40)+2)/40:nTjCellCut{1,2}(3,40)-(nTjCellCut{1,2}(3,40)+2)/40, 'Color',[79, 138, 205]./255,'LineStyle',':','LineWidth',1); 
plot3(ones(1,40)*nTjCellCut{1,3}(1,40), ones(1,40)*nTjCellCut{1,3}(2,40),-2:(nTjCellCut{1,3}(3,40)+2)/40:nTjCellCut{1,3}(3,40)-(nTjCellCut{1,3}(3,40)+2)/40, 'Color',[21, 67, 130]./255,'LineStyle',':','LineWidth',1); 

plot3(-2:(nTjCellCut{1,1}(1,40)+2)/40:nTjCellCut{1,1}(1,40)-(nTjCellCut{1,1}(1,40)+2)/40, ones(1,40)*nTjCellCut{1,1}(2,40),ones(1,40)*-2, 'Color',[39, 170, 225]./255,'LineStyle',':','LineWidth',1); 
plot3(-2:(nTjCellCut{1,2}(1,40)+2)/40:nTjCellCut{1,2}(1,40)-(nTjCellCut{1,2}(1,40)+2)/40, ones(1,40)*nTjCellCut{1,2}(2,40),ones(1,40)*-2, 'Color',[79, 138, 205]./255,'LineStyle',':','LineWidth',1); 
plot3(-2:(nTjCellCut{1,3}(1,40)+2)/40:nTjCellCut{1,3}(1,40)-(nTjCellCut{1,3}(1,40)+2)/40, ones(1,40)*nTjCellCut{1,3}(2,40),ones(1,40)*-2, 'Color',[21, 67, 130]./255,'LineStyle',':','LineWidth',1); 

plot3(ones(1,40)*nTjCellCut{1,1}(1,40), -1:(nTjCellCut{1,1}(2,40)+1)/40:nTjCellCut{1,1}(2,40)-(nTjCellCut{1,1}(2,40)+1)/40,ones(1,40)*-2, 'Color',[39, 170, 225]./255,'LineStyle',':','LineWidth',1); 
plot3(ones(1,40)*nTjCellCut{1,2}(1,40), -1:(nTjCellCut{1,2}(2,40)+1)/40:nTjCellCut{1,2}(2,40)-(nTjCellCut{1,2}(2,40)+1)/40,ones(1,40)*-2, 'Color',[79, 138, 205]./255,'LineStyle',':','LineWidth',1); 
plot3(ones(1,40)*nTjCellCut{1,3}(1,40), -1:(nTjCellCut{1,3}(2,40)+1)/40:nTjCellCut{1,3}(2,40)-(nTjCellCut{1,3}(2,40)+1)/40,ones(1,40)*-2, 'Color',[21, 67, 130]./255,'LineStyle',':','LineWidth',1); 

% plot3(ones(1,100)*valAtT(1,1),-2:3/100:1-3/100, ones(1,100)*-2,'Color',[39, 170, 225]./255,'LineStyle',':')
% plot3(ones(1,100)*valAtT(1,2),-2:3/100:1-3/100, ones(1,100)*-2,'Color',[79, 138, 205]./255,'LineStyle',':')
% plot3(ones(1,100)*valAtT(1,3),-2:3/100:1-3/100, ones(1,100)*-2,'Color',[21, 67, 130]./255,'LineStyle',':')
% 
% plot3(ones(1,100)*-2,-2:3.4/100:1.4-3.4/100, ones(1,100)*valAtT(3,1),'Color',[39, 170, 225]./255,'LineStyle',':')
% plot3(ones(1,100)*-2,-2:3.4/100:1.4-3.4/100, ones(1,100)*valAtT(3,2),'Color',[79, 138, 205]./255,'LineStyle',':')
% plot3(ones(1,100)*-2,-2:3.4/100:1.4-3.4/100, ones(1,100)*valAtT(3,3),'Color',[21, 67, 130]./255,'LineStyle',':')
% 
% plot3(ones(1,100)*4,-2:3.4/100:1.4-3.4/100, ones(1,100)*valAtT(3,1),'Color',[39, 170, 225]./255,'LineStyle',':')
% plot3(ones(1,100)*4,-2:3.4/100:1.4-3.4/100, ones(1,100)*valAtT(3,2),'Color',[79, 138, 205]./255,'LineStyle',':')
% plot3(ones(1,100)*4,-2:3.4/100:1.4-3.4/100, ones(1,100)*valAtT(3,3),'Color',[21, 67, 130]./255,'LineStyle',':')
% 
% plot3(-2:6/100:4-6/100,ones(1,100)*-1, ones(1,100)*valAtT(3,1),'Color',[39, 170, 225]./255,'LineStyle',':')
% plot3(-2:6/100:4-6/100,ones(1,100)*-1, ones(1,100)*valAtT(3,2),'Color',[79, 138, 205]./255,'LineStyle',':')
% plot3(-2:6/100:4-6/100,ones(1,100)*-1, ones(1,100)*valAtT(3,3),'Color',[21, 67, 130]./255,'LineStyle',':')
% 
% plot3(-2:6/100:4-6/100, ones(1,100)*valAtT(2,1), ones(1,100)*-2,'Color',[39, 170, 225]./255,'LineStyle',':')
% plot3(-2:6/100:4-6/100, ones(1,100)*valAtT(2,2), ones(1,100)*-2,'Color',[79, 138, 205]./255,'LineStyle',':')
% plot3(-2:6/100:4-6/100, ones(1,100)*valAtT(2,3), ones(1,100)*-2,'Color',[21, 67, 130]./255,'LineStyle',':')
% 
% plot3(-2:6/100:4-6/100, ones(1,100)*1.2, ones(1,100)*valAtT(3,1), 'Color',[39, 170, 225]./255,'LineStyle',':')
% plot3(-2:6/100:4-6/100, ones(1,100)*1.2, ones(1,100)*valAtT(3,2), 'Color',[79, 138, 205]./255,'LineStyle',':')
% plot3(-2:6/100:4-6/100, ones(1,100)*1.2, ones(1,100)*valAtT(3,3), 'Color',[21, 67, 130]./255,'LineStyle',':')

% try fill3 to insert planes corresponding to each reach amplitude tertile

xlim([-2 4]); ylim([-1 1.2]); zlim([-2 4])
yticks(-1:1:1)
print( fullfile(filePath,'Figure',strcat('projectionToETag-ITag-DistDimensions')), '-dpdf', '-bestfit', '-painters')

figure; % ITag-Dim x Distance-Dim Distance tertiles
hold on; 
% plot the max distance tertiles (trial-averaged)
plot2DneuralTrajAndEventMarkers( {nTjCell{1}([2,3],:)}, [39, 170, 225]./255, evtMarkers, evtCmap, 1, 5 ); 
plot2DneuralTrajAndEventMarkers( {nTjCell{2}([2,3],:)}, [79, 138, 205]./255, evtMarkers, evtCmap, 1, 5 ); 
plot2DneuralTrajAndEventMarkers( {nTjCell{3}([2,3],:)}, [21, 67, 130]./255, evtMarkers, evtCmap, 1, 5 ); 
xlim([-2 4]);
ylim([-2 4]);
hold off; 
print( fullfile(filePath,'Figure',strcat('projectionToITag-DistDimensions')), '-dpdf', '-bestfit')


%% X-Y plot neural Traj excursion amounts at max Js position along with ETag vs ITag dimensions 
nTjCellAllTrialsInhDim = cellfun(@(x) x(1,:), nTjCellAllTrialsInhMovDim, 'UniformOutput', false); 
nTjCellAllTrialsInhDim = squeeze(cell2mat(nTjCellAllTrialsInhDim(1,1,:)))'; % all trials projected onto the tagI dim across time bins
nTjCellAllTrialsInhDimNorm = nTjCellAllTrialsInhDim-mean(mean(nTjCellAllTrialsInhDim(:,1:10))); 

nTjCellAllTrialsExtDim = cellfun(@(x) x(1,:), nTjCellAllTrialsExtMovDim, 'UniformOutput', false); 
nTjCellAllTrialsExtDim = squeeze(cell2mat(nTjCellAllTrialsExtDim(1,1,:)))'; % all trials projected onto the tagI dim across time bins
nTjCellAllTrialsExtDimNorm = nTjCellAllTrialsExtDim-mean(mean(nTjCellAllTrialsExtDim(:,1:10))); 

for i = 1:sum(valTrialIdx)
    nTjCellAllTrialsInhDimNormMaxExc(i) = nTjCellAllTrialsInhDimNorm(i,maxReachPosIcrt(i)); 
end
clearvars i 

for i = 1:sum(valTrialIdx)
    nTjCellAllTrialsExtDimNormMaxExc(i) = nTjCellAllTrialsExtDimNorm(i,maxReachPosIcrt(i)); 
end
clearvars i 

% save tByt max excursion distance along with tagI, tagE dimensions with the discretized movement distance 
maxExcTagIEdim = [nTjCellAllTrialsInhDimNormMaxExc;  nTjCellAllTrialsExtDimNormMaxExc; maxP(:,3)'];  % save tByt max excursion distance along with tagI, tagE dimensions with the discretized movement distance 
save(fullfile(filePath,'IT01Ctx_121317_targetedDimReductionRegressionDistDimOnly'),'maxExcTagIEdim','ts','-append')

figure; 
hold on; 
plot(nTjCellAllTrialsInhDimNormMaxExc(maxExcTagIEdim(3,:)==1), nTjCellAllTrialsExtDimNormMaxExc(maxExcTagIEdim(3,:)==1), 'o', 'MarkerFaceColor', [39, 170, 225]./255, 'MarkerSize', 10, 'MarkerEdgeColor', 'none')
plot(nTjCellAllTrialsInhDimNormMaxExc(maxExcTagIEdim(3,:)==2), nTjCellAllTrialsExtDimNormMaxExc(maxExcTagIEdim(3,:)==2), 'o', 'MarkerFaceColor', [79, 138, 205]./255, 'MarkerSize', 10, 'MarkerEdgeColor', 'none')
plot(nTjCellAllTrialsInhDimNormMaxExc(maxExcTagIEdim(3,:)==3), nTjCellAllTrialsExtDimNormMaxExc(maxExcTagIEdim(3,:)==3), 'o', 'MarkerFaceColor', [21, 67, 130]./255, 'MarkerSize', 10, 'MarkerEdgeColor', 'none')
xlim([-4 10])
ylim([-4 10])
plot([-4 10], [-4 10],':k')
hold off; 






