function [nTrj] = visualizeNeuralTraj(filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, varargin)
%This function opens/visualizes the trial-by-trial neural
% population trajectories extracted using either pca or gpfa. 

cd(filePath)
p = parse_input_visualizeNeuralTraj( filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, varargin );
% p = parse_input_visualizeNeuralTraj( filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, {} );

%% load/organize neural trjectories pca- or gpfa-based
neuralTrajFile = dir(fullfile(filePath,fileNameNeuralTrj));
load(neuralTrajFile.name,'pcaResult')
nTrj.trjMat  = reshape([pcaResult.kern.seqTrain.xpost], pcaResult.p.Results.pcaDim, [], length(pcaResult.kern.seqTrain)); % neural trajectories in dim x timeBins x trials (e.g. 5x100x114)
nTrj.trialId = [pcaResult.kern.seqTrain.trialId]; % trial IDs for the neural trjectories

nTrj.relativeTimeBins = pcaResult.p.Results.timeRange(1):pcaResult.p.Results.pcaBinSize:pcaResult.p.Results.timeRange(2)-pcaResult.p.Results.pcaBinSize; % get timeBins relative to time 0 (e.g. reward) of the current PSTH
nTrj.neuralTrajFile = neuralTrajFile; % to keep track of the neural trajectory files
nTrj.pcaRez = pcaResult.kern; % pca info pc loadings, eigVals, expVar 

cd(fullfile(filePath,'Figure'))

%% plot rank-folded neural population trajectories
% plot principal components
%plot(pcaResult.kern.estParams.L(:,1)) % the 1st principal component
%plotEachDimVsTime(pcaResult.kern.seqTrain, 'xpost', pcaResult.binWidth);
trialByTimePCscore = squeeze(nTrj.trjMat(p.Results.dimPlotBy,:,:))'; % get the trial-by-timeBin PC score matrix for the current dimension
[~,sortTrI] = sortrows(max(abs(trialByTimePCscore),[],2),'descend'); % sort rows as a function of the projection score onto the dimension of interest

if strcmpi(p.Results.whatToPlot,'Trials')
    for t = 1:p.Results.numbTrials
        nTrjCell{t} = nTrj.trjMat(:,:,sortTrI(t));
    end
    clearvars t
elseif strcmpi(p.Results.whatToPlot,'Folds')
    trPerFold = floor(length(nTrj.trialId)/p.Results.numbFolds); % trial p.Results.trialFolds
    % get the nTraj mat sorted by the score as a result of projecting onto the dimension of interest
    for f = 0:p.Results.numbFolds-1
        if f<p.Results.numbFolds-1
            nTrjCell{f+1} = nanmean(nTrj.trjMat(:,:,sortTrI(f*trPerFold+1:(f+1)*trPerFold,1)),3); % take the trials of the current fold sorted by the max position
        elseif f==p.Results.numbFolds-1 % take all the remaining trials for the final fold
            nTrjCell{f+1} = nanmean(nTrj.trjMat(:,:,sortTrI(f*trPerFold+1:end,1)),3); % take the trials of the current fold sorted by the max position
        end
    end
    clearvars f
elseif strcmpi(p.Results.whatToPlot,'selectTrials')
     pcaRezTrials = [pcaResult.kern.seqTrain.trialId]; 
     for t = 1:length(p.Results.selectTrials) 
         if ismember(p.Results.selectTrials(t),pcaRezTrials)
             nTrjCell{t} = nTrj.trjMat(:,:,sortTrI(t));
         else
             warning(strcat(sprintf('Trial#%d',p.Results.selectTrials(t)),' was NOT found/plotted!'))
         end
     end
elseif strcmpi(p.Results.whatToPlot,'random')
    randTrials = randsample(length(sortTrI),p.Results.numbTrials);  % random sample trials
    for t = 1:p.Results.numbTrials
        nTrjCell{t} = nTrj.trjMat(:,:,randTrials(t));
    end
    clearvars t 
else
    error('Specify what/how to plot! Trials? random? Folds? selectTrials?')
end

eventMarkers = [1, arrayfun(@(x) find(x==nTrj.relativeTimeBins), eventMarkersRelative)]; % find the time bins corresponding to the event markers using the relative eventMarkers (e.g. [-1050 0]) given as the input

nTrjCmap = TNC_CreateRBColormapJP(length(nTrjCell), p.Results.trjCmap); % get the neural trajectory colormap
nTrjCmap = flip(nTrjCmap,1); 
eventMarkerCmap = cool(length(eventMarkers)+1); % get the colormap for all the events to be marked

%cellfun(@(x) nanmean(x(1:3,:,:),3), nTrj.sortedFoldMaxPos, 'UniformOutput', false); % take the trial-averaged neural trajectories of each fold based on the maxPos
plot3DneuralTrajAndEventMarkers( nTrjCell, nTrjCmap, eventMarkers, eventMarkerCmap, p.Results.lineWidth, p.Results.markerSize );  % plot the neural trajectories rank-ordered by the maxPos
title(strcat('trAvgSortedFoldMaxPos_',saveNameTag),'Interpreter', 'none')
if p.Results.saveNtrjFigs
    print(strcat(saveNameTag,'_trAvgSortedFoldMaxPos'),'-dpdf');
end

nTrj.nTrjCell = nTrjCell; 
nTrj.nTrjCmap = nTrjCmap; 

%% save
cd(filePath)
save(strcat(saveNameTag,'_visualize_nTrj'),'nTrj','p') % save the outcomes

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot3DneuralTrajAndEventMarkers( nTrjCell, neuralTrajCmap, eventMarkers, eventMarkerCmap, lineWidth, markerSize )
        figure;
        hold on;
        for fd = 1:length(nTrjCell)
            plot3(nTrjCell{fd}(1,:),nTrjCell{fd}(2,:),nTrjCell{fd}(3,:),'LineWidth',lineWidth,'color',neuralTrajCmap(fd,:))
            for evt = 1:length(eventMarkers)
                plot3(nTrjCell{fd}(1,eventMarkers(evt)),nTrjCell{fd}(2,eventMarkers(evt)),nTrjCell{fd}(3,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
            end
        end
        hold off;
        pbaspect([1 1 1]); grid on;
        %s=inputname(1); % take the input variable name as a string
        %title(strcat(s));
        xlabel('Dim1'); ylabel('Dim2'); zlabel('Dim3')
    end

    function p = parse_input_visualizeNeuralTraj( filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, vargs ) 
        %parse input, and extract name-value pairs for the main function 'corrNeuralTrajMovKinematics.m'
        
        default_whatToPlot = 'Trials'; 
        default_dimPlotBy = 1; % the default PC or GPFA dimension to sort the trials/folds by 
        default_numbTrials = 10; % # of trials to plot
        default_numbFolds = 4; % # of folds to plot
        default_selectTrials = 1:10; % default select trials to plot
        
        default_lineWidth  = 2;  % default lineWidth to be used to plot the neural trajectories
        default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
        default_saveNtrjFigs = false; % By default, there's no need to save the 3-d nTrj figures
        default_alpha = 0.05;  % default alpha value to just plot the significant correlations between neural population trajectories and the behavioral kinematics
        default_trjCmap = 'parula'; % the default trajectory colormap
        
        p = inputParser; % create parser object
        
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileNameNeuralTrj'); % fileName for neural population trajectories
        addRequired(p,'saveNameTag'); % saveName used to save the outcomes
        addRequired(p,'eventMarkersRelative'); % timing for specific events to be marked on each neural trajectory
        
        addParameter(p,'whatToPlot', default_whatToPlot) % three different plot modes: 'Trials', 'Folds', 'selectTrials'
        addParameter(p,'dimPlotBy', default_dimPlotBy) % the default PC or GPFA dimension to sort the trials/folds by
        addParameter(p,'numbTrials', default_numbTrials) % the # of trials to plot
        addParameter(p,'numbFolds', default_numbFolds) % the # of folds to plot
        addParameter(p,'selectTrials', default_selectTrials) % the select trials to plot
        
        addParameter(p,'lineWidth', default_lineWidth)   % the lineWidth to be used to plot nTrjs
        addParameter(p,'markerSize', default_markerSize) % the markerSzie to be used to plot nTrjs
        addParameter(p,'saveNtrjFigs', default_saveNtrjFigs) % By default, there's no need to save the 3-d nTrj figures
        addParameter(p,'alpha', default_alpha) % By default, the alpha level of 0.01 is used as the criterion for significant correlations
        addParameter(p,'trjCmap', default_trjCmap) % the default trajectory colormap 
        
        parse(p, filePath, fileNameNeuralTrj, saveNameTag, eventMarkersRelative, vargs{:})
        
    end

end

