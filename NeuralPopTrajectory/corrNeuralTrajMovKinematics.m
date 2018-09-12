function [nTrj, bTrj] = corrNeuralTrajMovKinematics(filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, varargin)
%corrNeuralTrajMovKinematics opens/loads the trial-by-trial neural
% population trajectories extracted using either pca or gpfa and movement kinematic
% variables saved in the 'behvariables.mat'. First, the correlation between the neural trajectory (e.g. PC scores across all dimensions)
% and the movement kinematics is examined  across all bins. Second, the trial-by-trial neural
% population trajectories are then divided into n folds (e.g. 4 folds)
% that are ranked based on a trial-by-trial movement kinematic variable - 1) max position, 2) max velocity, 3) lick counts.
% Finally, the neural population trajectories of all folds are drawn using the top three dims.

%filePath = '/Volumes/RAID2/parkj/NeuralData/ITphys/IT01_Ldms_M1_121317/Matfiles';
%fileNameNeuralTrj = 'IT01Str_121317_pca_reward_5D_50msBin.mat';
%fileNameBeh = 'BehVariables.mat';
%saveNameTag = 'IT01Str_121317_pca_reward_5D_50msBin';
%eventMarkersRelative = [-1050 0]; % time points to mark on the neural trajectories
%reachTimeWin = [-2500 0]; % timeWin to focus movement kinematics of
%lickTimeWin = [0 2000];   % timeWin to focus lick counts of
cd(filePath)
p = parse_input_corrNeuralTrajMovKinematics( filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, varargin );
% p = parse_input_corrNeuralTrajMovKinematics( filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, {} );

% for f = 1:p.Results.trialFolds
%     folds{1,f} = num2str(f); % get the fold labels in string
% end
% clearvars f

%% load/organize neural trjectories pca- or gpfa-based
neuralTrajFile = dir(fullfile(filePath,fileNameNeuralTrj));
load(neuralTrajFile.name,'pcaResult')
nTrj.trjMat = reshape([pcaResult.kern.seqTrain.xpost], pcaResult.p.Results.pcaDim, [], length(pcaResult.kern.seqTrain)); % neural trajectories in dim x timeBins x trials (e.g. 5x100x114)
nTrj.trialId = [pcaResult.kern.seqTrain.trialId]; % trial IDs for the neural trjectories
nTrj.relativeTimeBins = pcaResult.p.Results.timeRange(1):pcaResult.p.Results.pcaBinSize:pcaResult.p.Results.timeRange(2)-pcaResult.p.Results.pcaBinSize; % get timeBins relative to time 0 (e.g. reward) of the current PSTH
nTrj.reachTimeBins = nTrj.relativeTimeBins>=reachTimeWin(1) & nTrj.relativeTimeBins<reachTimeWin(2); % the time bins of interest for the behavioral kinematics (e.g. when reaches are most likely to occur) to narrow down the range of behavioral measures
nTrj.lickTimeBins  = nTrj.relativeTimeBins>=lickTimeWin(1) & nTrj.relativeTimeBins<lickTimeWin(2); % the time bins of interest for the behavioral kinematics (e.g. when reaches are most likely to occur) to narrow down the range of behavioral measures
nTrj.neuralTrajFile = neuralTrajFile; % to keep track of the neural trajectory files
nTrj.pcaRez = pcaResult.kern; % pca info pc loadings, eigVals, expVar 

%% load/organize behavioral data 'BehVariables.mat'
behFile = dir(fullfile(filePath,fileNameBeh));
load(behFile.name, 'ts', 'reach0','lick')

if ~contains(fileNameNeuralTrj,'reward') && contains(fileNameNeuralTrj,'reach')
    behTS = ts.reachStart(nTrj.trialId); % nTrj.trialId (*seqTrain.trialId) contains non-NaN trials only
elseif ~contains(fileNameNeuralTrj,'reach') && contains(fileNameNeuralTrj,'reward')
    behTS = ts.reward(nTrj.trialId); % nTrj.trialId (*seqTrain.trialId) contains non-NaN trials only
else
    error('Make sure if correct neural data has been input!')
end

lickBin = zeros(1,length(lick)); % licks
lickBin(ts.lick) = 1; % licks

valTrialIdx = zeros(length(behTS),1);
for t = 1:length(behTS) % increment trials, take the trial-by-trial position/velocity data and bin them to match the neural trjectories (e.g. 50 ms)
    if max(behTS(t)+pcaResult.p.Results.timeRange)<=length(reach0) && min(behTS(t)+pcaResult.p.Results.timeRange)>0
        valTrialIdx(t)=1;
        timeWin = behTS(t)+pcaResult.p.Results.timeRange; % the time window, e.g. -3 to 2 sec relative to the behavioral timestamp
        bTrj(t).reachPos = smooth(decimate(reach0(timeWin(1):timeWin(2)-1),pcaResult.p.Results.pcaBinSize),3)'; % get decimated behavioral trjectories on the same timescale of the neural trjectories
        bTrj(t).reachVel = smooth(diff([bTrj(t).reachPos(1) bTrj(t).reachPos]),3)'; % get reach velocities (from reach0)
        bTrj(t).maxReachPos = max(bTrj(t).reachPos(nTrj.reachTimeBins)); % max reach position within the time window of interest
        bTrj(t).maxReachVel = max(bTrj(t).reachVel(nTrj.reachTimeBins)); % max reach velocity within the time window of interest
        bTrj(t).trialId = t;
        bTrj(t).lick = bin1msSpkCountMat(lickBin(timeWin(1):timeWin(2)-1),pcaResult.p.Results.pcaBinSize,pcaResult.p.Results.pcaBinSize); %
        bTrj(t).lickCount =  sum(bTrj(t).lick(nTrj.lickTimeBins)); % lick Counts within the lickTimeBins
    else
    end
end
clearvars t
%imagesc((reshape([bTrj.lick]', length(nTrj.relativeTimeBins), length(nTrj.trialId)))');

%% correlation between neural population trjectories and movement kinematics
for dim = 1:size(nTrj.trjMat,1) % increment dimensions
    trialByTimePCscore =  squeeze(nTrj.trjMat(dim,:,valTrialIdx==1))'; % get the trial-by-timeBin PC score matrix for the current dimension
    % correlation movement kinematics & neural trjectories
    [nTrj.rMaxPos(:,dim),nTrj.pMaxPos(:,dim)] = corr(trialByTimePCscore,[bTrj.maxReachPos]'); % correlation between PC score in the current dim across all time bins and maxReachPos
    [nTrj.rMaxVel(:,dim),nTrj.pMaxVel(:,dim)] = corr(trialByTimePCscore,[bTrj.maxReachVel]'); % correlation between PC score in the current dim across all time bins and maxReachVel
    [rPosBbyBtotal,pPosBbyBtotal] = corr(trialByTimePCscore,reshape([bTrj.reachPos]',[],sum(valTrialIdx))');  % correlation between PC score in the current dim across all time bins and timeBin-by-timeBin reachPos
    nTrj.rPosBbyB(:,dim) = diag(rPosBbyBtotal); % take rho of the matching time bins
    nTrj.pPosBbyB(:,dim) = diag(pPosBbyBtotal); % take pVal of the matching time bins
    [rVelBbyBtotal,pVelBbyBtotal] = corr(trialByTimePCscore,reshape([bTrj.reachVel]',[],sum(valTrialIdx))');  % corr timeBin-by-timeBin
    nTrj.rVelBbyB(:,dim) = diag(rVelBbyBtotal); % take rho of the matching time bins
    nTrj.pVelBbyB(:,dim) = diag(pVelBbyBtotal); % take pVal of the matching time bins
    % correlation lick counts & neural trjectories
    [nTrj.rLick(:,dim),nTrj.pLick(:,dim)] = corr(trialByTimePCscore,[bTrj.lickCount]'); % correlation between the current
end
clearvars dim

cd(fullfile(filePath,'Figure'))

cmap = TNC_CreateRBColormap(100,'rb'); % generate a colormap for imagesc psth
figure;
rMaxPosSig = zeros(size(nTrj.pMaxPos,1), size(nTrj.pMaxPos,2));
rMaxPosSig(nTrj.pMaxPos<p.Results.alpha) = nTrj.rMaxPos(nTrj.pMaxPos<p.Results.alpha);
imagescJP(rMaxPosSig',cmap,[-0.8 0.8]);
pbaspect([1 1 1]); colorbar;
title(strcat('CorrPCscrDimMaxPos_',saveNameTag),'Interpreter', 'none');
set(gca,'yTick',1:1:size(nTrj.rMaxPos',1)); % plot the corr between nTrj and maxPos dim-by-dim
print(strcat(saveNameTag,'_corrPCscrDimMaxPos'),'-dpdf');

figure;
rMaxVelSig = zeros(size(nTrj.pMaxVel,1), size(nTrj.pMaxVel,2));
rMaxVelSig(nTrj.pMaxVel<p.Results.alpha) = nTrj.rMaxVel(nTrj.pMaxVel<p.Results.alpha);
imagescJP(rMaxVelSig',cmap,[-0.8 0.8]);
pbaspect([1 1 1]); colorbar;
title(strcat('CorrPCscrDimMaxVel_',saveNameTag),'Interpreter', 'none');
set(gca,'yTick',1:1:size(nTrj.rMaxVel',1)); % plot the corr between nTrj and maxVel dim-by-dim
print(strcat(saveNameTag,'_corrPCscrDimMaxVel'),'-dpdf');

figure;
rLickSig = zeros(size(nTrj.pLick,1), size(nTrj.pLick,2));
rLickSig(nTrj.pLick<p.Results.alpha) = nTrj.rLick(nTrj.pLick<p.Results.alpha);
imagescJP(rLickSig',cmap,[-0.8 0.8]);
pbaspect([1 1 1]); colorbar;
title(strcat('CorrPCscrDimLickCnt_',saveNameTag),'Interpreter', 'none');
set(gca,'yTick',1:1:size(nTrj.rLick',1)); % plot the corr between nTrj and lickCount dim-by-dim
print(strcat(saveNameTag,'_corrPCscrDimLickCnt'),'-dpdf');

%% plot rank-folded neural population trajectories
% plot principal components
%plot(pcaResult.kern.estParams.L(:,1)) % the 1st principal component
%plotEachDimVsTime(pcaResult.kern.seqTrain, 'xpost', pcaResult.binWidth);

% plot trajectories
trFolds = floor(length(bTrj)/p.Results.trialFolds); % trial p.Results.trialFolds
sortMaxPos = sortrows([[bTrj.maxReachPos];[bTrj.trialId]]',-1); % sorted max reach positions in a descending order
sortMaxVel = sortrows([[bTrj.maxReachVel];[bTrj.trialId]]',-1); % sorted max reach velocities
sortLickCnt = sortrows([[bTrj.lickCount];[bTrj.trialId]]',-1);  % sorted trial lick counts

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 0:p.Results.trialFolds-1
    if f<p.Results.trialFolds-1
        nTrj.sortedFoldMaxPos{f+1} = nTrj.trjMat(:,:,sortMaxPos(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by the max position
        nTrj.sortedFoldMaxVel{f+1} = nTrj.trjMat(:,:,sortMaxVel(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by the max velocity
        nTrj.sortedFoldLickCnt{f+1} = nTrj.trjMat(:,:,sortLickCnt(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by the lick counts
    elseif f==p.Results.trialFolds-1 % take all the remaining trials for the final fold
        nTrj.sortedFoldMaxPos{f+1} = nTrj.trjMat(:,:,sortMaxPos(f*trFolds+1:end,2)); % take the trials of the current fold sorted by the max position
        nTrj.sortedFoldMaxVel{f+1} = nTrj.trjMat(:,:,sortMaxVel(f*trFolds+1:end,2)); % take the trials of the current fold sorted by the max velocity
        nTrj.sortedFoldLickCnt{f+1} = nTrj.trjMat(:,:,sortLickCnt(f*trFolds+1:end,2)); % take the trials of the current fold sorted by the lick counts
    end
end
clearvars f

eventMarkers = [1, arrayfun(@(x) find(x==nTrj.relativeTimeBins), eventMarkersRelative)]; % find the time bins corresponding to the event markers using the relative eventMarkers (e.g. [-1050 0]) given as the input
neuralTrajCmap  = summer(p.Results.trialFolds); % get the colormap for fold-by-fold trial-averaged neural trajectories
eventMarkerCmap = cool(length(eventMarkers)+1); % get the colormap for all the events to be marked

trAvgSortedFoldMaxPos = cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldMaxPos, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos
trAvgSortedFoldMaxVel = cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldMaxVel, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos
trAvgSortedFoldLickCnt = cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldLickCnt, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos

plot3DneuralTrajAndEventMarkers( trAvgSortedFoldMaxPos, neuralTrajCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize );  % plot the neural trajectories rank-ordered by the maxPos
title(strcat('trAvgSortedFoldMaxPos_',saveNameTag),'Interpreter', 'none')
if p.Results.saveNtrjFigs
    print(strcat(saveNameTag,'_trAvgSortedFoldMaxPos'),'-dpdf');
end

plot3DneuralTrajAndEventMarkers( trAvgSortedFoldMaxVel, neuralTrajCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize );  % plot the neural trajectories rank-ordered by the maxVel
title(strcat('trAvgSortedFoldMaxVel_',saveNameTag),'Interpreter', 'none')
if p.Results.saveNtrjFigs
    print(strcat(saveNameTag,'_trAvgSortedFoldMaxVel'),'-dpdf');
end

plot3DneuralTrajAndEventMarkers( trAvgSortedFoldLickCnt, neuralTrajCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize ); % plot the neural trajectories rank-ordered by the lickCnt
title(strcat('trAvgSortedFoldLickCnt_',saveNameTag),'Interpreter', 'none')
if p.Results.saveNtrjFigs
    print(strcat(saveNameTag,'_trAvgSortedFoldLickCnt'),'-dpdf');
end

%% save
cd(filePath)
save(strcat(saveNameTag,'_nTrj_bTrj_Corr'),'nTrj','bTrj','p') % save the outcomes

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot3DneuralTrajAndEventMarkers( neuralTrajCell, neuralTrajCmap, eventMarkers, eventMarkerCmap, lineWidth, markerSize )
        figure;
        hold on;
        for fd = 1:length(neuralTrajCell)
            plot3(neuralTrajCell{fd}(1,:),neuralTrajCell{fd}(2,:),neuralTrajCell{fd}(3,:),'LineWidth',lineWidth,'color',neuralTrajCmap(fd,:))
            for evt = 1:length(eventMarkers)
                plot3(neuralTrajCell{fd}(1,eventMarkers(evt)),neuralTrajCell{fd}(2,eventMarkers(evt)),neuralTrajCell{fd}(3,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
                plot3(neuralTrajCell{fd}(1,eventMarkers(evt)),neuralTrajCell{fd}(2,eventMarkers(evt)),neuralTrajCell{fd}(3,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
                plot3(neuralTrajCell{fd}(1,eventMarkers(evt)),neuralTrajCell{fd}(2,eventMarkers(evt)),neuralTrajCell{fd}(3,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
            end
        end
        hold off;
        pbaspect([1 1 1]); grid on;
        %s=inputname(1); % take the input variable name as a string
        %title(strcat(s));
        xlabel('Dim1'); ylabel('Dim2'); zlabel('Dim3')
    end



    function p = parse_input_corrNeuralTrajMovKinematics( filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        %parse input, and extract name-value pairs for the main function 'corrNeuralTrajMovKinematics.m'
        
        default_trialFolds = 4;  % # of p.Results.trialFolds to classify the total trials by the rank of a behavioral variable
        default_lineWidth  = 2;  % default lineWidth to be used to plot the neural trajectories
        default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
        default_saveNTjcFigs = false; % By default, there's no need to save the 3-d nTrj figures
        default_alpha = 0.05; % default alpha value to just plot the significant correlations between neural population trajectories and the behavioral kinematics
        
        p = inputParser; % create parser object
        
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileNameNeuralTrj'); % fileName for neural population trajectories
        addRequired(p,'fileNameBeh'); % fileName for movement kinematic variables
        addRequired(p,'saveNameTag'); % saveName used to save the outcomes
        addRequired(p,'eventMarkersRelative'); % timing for specific events to be marked on each neural trajectory
        addRequired(p,'reachTimeWin'); % timeWin in which the movement kinematics are concerned
        addRequired(p,'lickTimeWin');  % timeWin in which lick counts are concerned
        
        addParameter(p,'trialFolds', default_trialFolds) % the # of p.Results.trialFolds to divide the whole trials into
        addParameter(p,'lineWidth', default_lineWidth)   % the lineWidth to be used to plot nTrjs
        addParameter(p,'markerSize', default_markerSize) % the markerSzie to be used to plot nTrjs
        addParameter(p,'saveNtrjFigs', default_saveNTjcFigs) % By default, there's no need to save the 3-d nTrj figures
        addParameter(p,'alpha', default_alpha) % By default, the alpha level of 0.01 is used as the criterion for significant correlations
        
        parse(p, filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, vargs{:})
        
    end

end

