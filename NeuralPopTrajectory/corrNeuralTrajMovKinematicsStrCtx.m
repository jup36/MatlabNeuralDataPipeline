function corrNeuralTrajMovKinematicsStrCtx(filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, varargin)
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

%% load/organize neural trjectories pca- or gpfa-based
neuralTrajFile = dir(fullfile(filePath,fileNameNeuralTrj));
load(neuralTrajFile.name,'pcaResult')
nTrj.trjMat = reshape([pcaResult.kern.seqTrain.xpost], pcaResult.p.Results.pcaDim, [], length(pcaResult.kern.seqTrain)); % neural trajectories in dim x timeBins x trials (e.g. 5x100x114)
nTrj.trialId = [pcaResult.kern.seqTrain.trialId]; % trial IDs for the neural trjectories
nTrj.relativeTimeBins = pcaResult.p.Results.timeRange(1):pcaResult.p.Results.pcaBinSize:pcaResult.p.Results.timeRange(2)-pcaResult.p.Results.pcaBinSize; % get timeBins relative to time 0 (e.g. reward) of the current PSTH
nTrj.reachTimeBins = nTrj.relativeTimeBins>=reachTimeWin(1) & nTrj.relativeTimeBins<reachTimeWin(2); % the time bins of interest for the behavioral kinematics (e.g. when reaches are most likely to occur) to narrow down the range of behavioral measures
nTrj.lickTimeBins  = nTrj.relativeTimeBins>=lickTimeWin(1) & nTrj.relativeTimeBins<lickTimeWin(2); % the time bins of interest for the behavioral kinematics (e.g. when reaches are most likely to occur) to narrow down the range of behavioral measures
nTrj.neuralTrajFile = neuralTrajFile; % to keep track of the neural trajectory files

%% get separate neural trajectories for each popultion (STR & CTX) using the nTrj.trjMat
prjHiDimNtrj  = ([pcaResult.kern.seqTrain.xpost]'*pcaResult.kern.estParams.L')'; % project back to the high dimensional space

strIdx = cell2mat(pcaResult.isStriatum); % index for striatal units
ctxIdx = ~strIdx; % index for cortical units

% project neural trajectories of each region back to the low-dim space
strPrjNtrj = (prjHiDimNtrj(strIdx(pcaResult.unitIdx),:)'*pcaResult.kern.estParams.L(strIdx(pcaResult.unitIdx),:))'; % project the striatal population activity separately back to the low-dim space
ctxPrjNtrj = (prjHiDimNtrj(ctxIdx(pcaResult.unitIdx),:)'*pcaResult.kern.estParams.L(ctxIdx(pcaResult.unitIdx),:))'; % project the cortical population activity separately back to the low-dim space

nTrj.prjTrjMat{1} = reshape(strPrjNtrj, pcaResult.p.Results.pcaDim, [], length(pcaResult.kern.seqTrain)); % reshape to unit x timeBin x trial mat
nTrj.prjTrjMat{2} = reshape(ctxPrjNtrj, pcaResult.p.Results.pcaDim, [], length(pcaResult.kern.seqTrain)); % reshape to unit x timeBin x trial mat

%% load/organize behavioral data 'BehVariables.mat'
behFile = dir(fullfile(filePath,fileNameBeh));
load(behFile.name, 'ts', 'reach0','lick')
behTS = ts.reward(nTrj.trialId); % nTrj.trialId (*seqTrain.trialId) contains non-NaN trials only
lickBin = zeros(1,length(lick)); % licks
lickBin(ts.lick) = 1; % licks

for t = 1:length(behTS) % increment trials, take the trial-by-trial position/velocity data and bin them to match the neural trjectories (e.g. 50 ms)
    timeWin = behTS(t)+pcaResult.p.Results.timeRange; % the time window, e.g. -3 to 2 sec relative to the behavioral timestamp
    bTrj(t).reachPos = smooth(decimate(reach0(timeWin(1):timeWin(2)-1),pcaResult.p.Results.pcaBinSize),3)'; % get decimated behavioral trjectories on the same timescale of the neural trjectories
    bTrj(t).reachVel = smooth(diff([bTrj(t).reachPos(1) bTrj(t).reachPos]),3)'; % get reach velocities (from reach0)
    bTrj(t).maxReachPos = max(bTrj(t).reachPos(nTrj.reachTimeBins)); % max reach position within the time window of interest
    bTrj(t).maxReachVel = max(bTrj(t).reachVel(nTrj.reachTimeBins)); % max reach velocity within the time window of interest
    bTrj(t).trialId = t;
    bTrj(t).lick = bin1msSpkCountMat(lickBin(timeWin(1):timeWin(2)-1),pcaResult.p.Results.pcaBinSize,pcaResult.p.Results.pcaBinSize); %
    bTrj(t).lickCount =  sum(bTrj(t).lick(nTrj.lickTimeBins)); % lick Counts within the lickTimeBins
end
clearvars t
%imagesc((reshape([bTrj.lick]', length(nTrj.relativeTimeBins), length(nTrj.trialId)))');

%% correlation between neural population trjectories and movement kinematics
cd(fullfile(filePath,'Figure'))
cmap = TNC_CreateRBColormap(100,'rb'); % generate a colormap to plot rho
for pop = 1:length(nTrj.prjTrjMat) % increment neural populations (e.g. str and ctx)
    for dim = 1:size(nTrj.prjTrjMat{pop},1) % increment dimensions
        trialByTimePCscore =  squeeze(nTrj.prjTrjMat{pop}(dim,:,:))'; % get the trial-by-timeBin PC score matrix for the current dimension
        % correlation movement kinematics & neural trjectories
        [nTrj.rMaxPos{pop}(:,dim),nTrj.pMaxPos{pop}(:,dim)] = corr(trialByTimePCscore,[bTrj.maxReachPos]'); % correlation between PC score in the current dim across all time bins and maxReachPos
        [nTrj.rMaxVel{pop}(:,dim),nTrj.pMaxVel{pop}(:,dim)] = corr(trialByTimePCscore,[bTrj.maxReachVel]'); % correlation between PC score in the current dim across all time bins and maxReachVel
        [rPosBbyBtotal,pPosBbyBtotal] = corr(trialByTimePCscore,reshape([bTrj.reachPos]',[],length(behTS))');  % correlation between PC score in the current dim across all time bins and timeBin-by-timeBin reachPos
        nTrj.rPosBbyB{pop}(:,dim) = diag(rPosBbyBtotal); % take rho of the matching time bins
        nTrj.pPosBbyB{pop}(:,dim) = diag(pPosBbyBtotal); % take pVal of the matching time bins
        [rVelBbyBtotal,pVelBbyBtotal] = corr(trialByTimePCscore,reshape([bTrj.reachVel]',[],length(behTS))');  % corr timeBin-by-timeBin
        nTrj.rVelBbyB{pop}(:,dim) = diag(rVelBbyBtotal); % take rho of the matching time bins
        nTrj.pVelBbyB{pop}(:,dim) = diag(pVelBbyBtotal); % take pVal of the matching time bins
        % correlation lick counts & neural trjectories
        [nTrj.rLick{pop}(:,dim),nTrj.pLick{pop}(:,dim)] = corr(trialByTimePCscore,[bTrj.lickCount]'); % correlation between the current
    end
    
    % image rMaxPos dim-by-dim across time bins and save it
    figure;
    imagescJP(nTrj.rMaxPos{pop}',cmap,[-0.8 0.8]); % image the rMaxPos
    pbaspect([1 1 1]);
    figTtlFmt = 'corrPCscrDimMaxPosNeuralPop#%d'; % figure title format
    figTtl = sprintf(figTtlFmt,pop); % format figure title into string
    title(figTtl); % figure title
    colorbar;
    set(gca,'yTick',1:1:size(nTrj.rMaxPos{pop}',1)); % plot the corr between nTrj and maxPos dim-by-dim
    print(strcat(saveNameTag,'_',figTtl),'-dpdf');
    
    % image rMaxVel dim-by-dim across time bins and save it
    figure;
    imagescJP(nTrj.rMaxVel{pop}',cmap,[-0.8 0.8]); % image the rMaxPos
    pbaspect([1 1 1]);
    figTtlFmt = 'corrPCscrDimMaxVelNeuralPop#%d'; % figure title format
    figTtl = sprintf(figTtlFmt,pop); % format figure title into string
    title(figTtl); % figure title
    colorbar;
    set(gca,'yTick',1:1:size(nTrj.rMaxVel{pop}',1)); % plot the corr between nTrj and maxPos dim-by-dim
    print(strcat(saveNameTag,'_',figTtl),'-dpdf');
    
    % image rLick dim-by-dim across time bins and save it
    figure;
    imagescJP(nTrj.rLick{pop}',cmap,[-0.8 0.8]); % image the rMaxPos
    pbaspect([1 1 1]);
    figTtlFmt = 'corrPCscrDimLickCntNeuralPop#%d'; % figure title format
    figTtl = sprintf(figTtlFmt,pop); % format figure title into string
    title(figTtl); % figure title
    colorbar;
    set(gca,'yTick',1:1:size(nTrj.rLick{pop}',1)); % plot the corr between nTrj and maxPos dim-by-dim
    print(strcat(saveNameTag,'_',figTtl),'-dpdf');
end
clearvars dim pop

%% plot rank-folded neural population trajectories
% plot principal components
%plot(pcaResult.kern.estParams.L(:,1)) % the 1st principal component
%plotEachDimVsTime(pcaResult.kern.seqTrain, 'xpost', pcaResult.binWidth);

% plot trajectories rank-ordered by a behavioral variable
trFolds = floor(length(bTrj)/p.Results.trialFolds); % trial p.Results.trialFolds
sortMaxPos = sortrows([[bTrj.maxReachPos];[bTrj.trialId]]',-1); % sorted max reach positions in a descending order
sortMaxVel = sortrows([[bTrj.maxReachVel];[bTrj.trialId]]',-1); % sorted max reach velocities
sortLickCnt = sortrows([[bTrj.lickCount];[bTrj.trialId]]',-1);  % sorted trial lick counts

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for pop = 1:length(nTrj.prjTrjMat) % increment neural populations (e.g. str and ctx)
    for f = 0:p.Results.trialFolds-1 % increment folds
        if f<p.Results.trialFolds-1
            nTrj.sortedFoldMaxPos{pop}{f+1} = nTrj.prjTrjMat{pop}(:,:,sortMaxPos(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by the max position
            nTrj.sortedFoldMaxVel{pop}{f+1} = nTrj.prjTrjMat{pop}(:,:,sortMaxVel(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by the max velocity
            nTrj.sortedFoldLickCnt{pop}{f+1} = nTrj.prjTrjMat{pop}(:,:,sortLickCnt(f*trFolds+1:(f+1)*trFolds,2)); % take the trials of the current fold sorted by the lick counts
        elseif f==p.Results.trialFolds-1 % take all the remaining trials for the final fold
            nTrj.sortedFoldMaxPos{pop}{f+1} = nTrj.prjTrjMat{pop}(:,:,sortMaxPos(f*trFolds+1:end,2)); % take the trials of the current fold sorted by the max position
            nTrj.sortedFoldMaxVel{pop}{f+1} = nTrj.prjTrjMat{pop}(:,:,sortMaxVel(f*trFolds+1:end,2)); % take the trials of the current fold sorted by the max velocity
            nTrj.sortedFoldLickCnt{pop}{f+1} = nTrj.prjTrjMat{pop}(:,:,sortLickCnt(f*trFolds+1:end,2)); % take the trials of the current fold sorted by the lick counts
        end
    end
end
clearvars f pop

eventMarkers = [1, arrayfun(@(x) find(x==nTrj.relativeTimeBins), eventMarkersRelative)]; % find the time bins corresponding to the event markers using the relative eventMarkers (e.g. [-1050 0]) given as the input
neuralTrajCmap  = summer(p.Results.trialFolds); % get the colormap for fold-by-fold trial-averaged neural trajectories
eventMarkerCmap = cool(length(eventMarkers)+1); % get the colormap for all the events to be marked

for pop = 1:length(nTrj.prjTrjMat) % increment neural populations (e.g. str and ctx)
    nTrj.trAvgSortedFoldMaxPos{pop} = cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldMaxPos{pop}, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos
    nTrj.trAvgSortedFoldMaxVel{pop} = cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldMaxVel{pop}, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos
    nTrj.trAvgSortedFoldLickCnt{pop} = cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldLickCnt{pop}, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos
    
    % plot trajectories rank-ordered by maxPos
    plot3DneuralTrajAndEventMarkers( nTrj.trAvgSortedFoldMaxPos{pop}, neuralTrajCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize );  % plot the neural trajectories rank-ordered by the maxPos
    figTtlFmt1 = 'trAvgSortedFoldMaxPosNeuralPop#%d'; % figure title format
    figTtl1 = sprintf(figTtlFmt1,pop); % format figure title into string
    title(figTtl1); % figure title
    if p.Results.saveNtrjFigs
        print(strcat(saveNameTag,'_',figTtl1),'-dpdf');
    end
    
    % plot trajectories rank-ordered by maxVel
    plot3DneuralTrajAndEventMarkers( nTrj.trAvgSortedFoldMaxVel{pop}, neuralTrajCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize );  % plot the neural trajectories rank-ordered by the maxVel
    figTtlFmt2 = 'trAvgSortedFoldMaxVelNeuralPop#%d'; % figure title format
    figTtl2 = sprintf(figTtlFmt2,pop); % format figure title into string
    title(figTtl2); % figure title
    if p.Results.saveNtrjFigs
        print(strcat(saveNameTag,'_trAvgSortedFoldMaxVel'),'-dpdf');
    end
    
    % plot trajectories rank-ordered by lickCnt
    plot3DneuralTrajAndEventMarkers( nTrj.trAvgSortedFoldLickCnt{pop}, neuralTrajCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize ); % plot the neural trajectories rank-ordered by the lickCnt
    figTtlFmt3 = 'trAvgSortedFoldLickCntNeuralPop#%d'; % figure title format
    figTtl3 = sprintf(figTtlFmt3,pop); % format figure title into string
    title(figTtl3); % figure title
    if p.Results.saveNtrjFigs
        print(strcat(saveNameTag,'_trAvgSortedFoldLickCnt'),'-dpdf');
    end
end
clearvars pop

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
        %title(s);
        xlabel('Dim1'); ylabel('Dim2'); zlabel('Dim3')
    end

    function p = parse_input_corrNeuralTrajMovKinematics( filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
        %parse input, and extract name-value pairs for the main function 'corrNeuralTrajMovKinematics.m'
        
        default_trialFolds = 4;  % # of p.Results.trialFolds to classify the total trials by the rank of a behavioral variable
        default_lineWidth  = 2;  % default lineWidth to be used to plot the neural trajectories
        default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
        default_saveNTjcFigs = false; % By default, there's no need to save the 3-d nTrj figures
        
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
        
        parse(p, filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, vargs{:})
        
    end

end

