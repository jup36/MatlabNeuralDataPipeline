function [ pcaResult, pcaDat ] = runPCAtrialAveragedTrajectories( filePath, fileName, saveNameTag, varargin )
%runPCA prepares neural population data in a format to extract low-dim trial-by-trial neural population trajectories 
% using Byron's neuralTraj package. 
% The input structure should provide trial x time spikecount matrix (binned in 1-ms). 
% The output structure 'pcaDat' must be in the dimension of Trials with 2 fields (trialId, spikes). 
% Each entry (trial) of the field 'spikes' must be Neuron-by-Timebins matrix. 
% pc loadings (weights) saved as 'pcaResult.kern.estParams.L'
% Important to note that PCA here extracts principal components of the
% neuronal activity covariance structure. The feature for PCA is different
% neurons in the population, of which binned spike counts observed across
% all time bins and trials (observation). PCA extract principal components
% of the neuronal activity covariance structure. 

p = parse_input_runPCA( filePath, fileName, saveNameTag, varargin );
% p = parse_input_runPCA( filePath, fileName, saveNameTag, {} ); % use this when running line-by-line

cd(p.Results.filePath)

pcaFileFolder = fullfile(p.Results.filePath,'mat_results','run*'); % file folder where the pca results are saved  
pcaRunNumber  = length(dir(pcaFileFolder))+1; % get the right pcaRunNumber to run neuralTraj.m, and save the result in the 'mat_results" folder under the current filePath. 

load(p.Results.fileName, 'S', 'pc'); % load the raw data structure, and results of the pca analysis using 'pcaPSTH.m' 
%load('BehVariables.mat','pos1')

% Check the existence of files
if exist('S','var') && ~exist('pc','var') 
    error('The variable pc is missing!!')
elseif ~exist('S','var') && exist('pc','var') 
    error('The variable S is missing!!')
elseif ~exist('S','var') && ~exist('pc','var') 
    error('Both S & pc are missing!!')
end

% Sanity check for the number of trials for each unit - check all units' spkCountMat have the same number of trials
if length(unique(cellfun(@length, S(:).SpkCountMat))) == 1 
   numbTrial = unique(cellfun(@length,S(:).SpkCountMat));
else
   error('There are units with different number of trials!!!')
end

valTimeBins = find( S.params{1}.binEdges1ms>=p.Results.timeRange(1) & S.params{1}.binEdges1ms<p.Results.timeRange(2) ); % find bins within the time range
unitMeanFR  = zeros(size(S.SpkCountMat,1),1); % unit mean FR rates within the valid time range to exclude low and high FR (e.g. < 2Hz & >50Hz) units
unitTimeTrial = zeros(size(S.SpkCountMat,1), length(valTimeBins), numbTrial); % unit x time(1ms bin) x trial matrix 

% organize the unitTime mat (trial-average)
for u = 1:size(S.SpkCountMat,1) % increment units
    tmpUnitTrialTimeMat  = full(cell2mat(S.SpkCountMat{u})); % trial x time spike count mat
    unitTimeTrial(u,:,:) = permute(tmpUnitTrialTimeMat(:,valTimeBins),[3 2 1]); % permute to get the unitTimeTrial mat
    unitMeanFR(u) = sum(sum(unitTimeTrial(u,:,:)))/(numbTrial*length(valTimeBins))*1000; % get the mean FR(Hz) of the unit
end
clearvars u 

% get rid of one unit of the pairs with high cross-correlation (remove duplicate/split units)
unitTimeAcrossTrials = reshape(unitTimeTrial, [size(unitTimeTrial,1), size(unitTimeTrial,2)*size(unitTimeTrial,3)]); % reshape the unitTimeTrial matrix to get unitTimeAcrossTrials (concatenate across trials)
unitTimeTrialAverage = nanmean(unitTimeTrial,3); % take average across trials
% determine how to block trials here based on which movement kinematic feature
% S.nonNaNtrialId contains 

xcorUnitIdx = getUnitOfHighXcorr( unitTimeAcrossTrials, p.Results.xcorThresholdPer );    
unitIdxPCA = unitMeanFR > p.Results.frHighPass & unitMeanFR < p.Results.frLowPass & xcorUnitIdx; % units to be used for pca
. 


if isfield(S,'nonNaNtrialId') % in case the nanNaNtrialId has been already sorted by the previous function 'pcaPSTHtwoPopulations.m'
    for tr = 1:length(S.nonNaNtrialId) % increment trials
        pcaDat(tr).trialId = S.nonNaNtrialId(tr,1); % trial ID
        pcaDat(tr).spikes = unitTimeTrial(unitIdxPCA,:,S.nonNaNtrialId(tr,1)); % neuron x timeBin mat for the trial, to include pca-sorted units only, pc.pcUnitSort(pc.pcUnitSort(:,3)==1,1)
    end
else
    for tr = 1:numbTrial % increment trials
        pcaDat(tr).trialId = tr; % trial ID
        pcaDat(tr).spikes = unitTimeTrial(unitIdxPCA,:,tr); % neuron x timeBin mat for the trial, to include pca-sorted units only, pc.pcUnitSort(pc.pcUnitSort(:,3)==1,1)
    end
    
end
clearvars tr

%% Extract neural trajectories 
pcaResult = neuralTraj(pcaRunNumber, pcaDat, 'datFormat', 'spikes', 'method', 'pca', 'xDim', p.Results.pcaDim, 'binWidth', p.Results.pcaBinSize, 'kernSDList', p.Results.pcaKernSD);

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.

% NOTES:
% - These function calls are computationally demanding.  Cross-validation 
%   takes a long time because a separate model has to be fit for each 
%   state dimensionality and each cross-validation fold.

if p.Results.fullCrossVal % cross-validation is not run by default
    
    for crossValXdim = [3 5 8]
        %neuralTraj(pcaRunNumber+1, pcaDat, 'method',  'pca', 'xDim', crossValXdim, 'numFolds', p.Results.crossValFolds);
        neuralTraj(pcaRunNumber+1, pcaDat, 'method', 'pca', 'xDim', crossValXdim, 'numFolds', p.Results.crossValFolds);
    end
    
    % Plot prediction error versus state dimensionality or versus kernSD.
        % Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
%     crossValKernSD = 50; % select kernSD for two-stage methods
%     plotPredErrorVsDim(pcaRunNumber+1, crossValKernSD);
%     crossValXDim = 5; % select state dimensionality
%     plotPredErrorVsKernSD(pcaRunNumber+1, crossValXDim);
    
end

%% Save pcaResult separately (in addition to the automatic save by 'neuralTraj')
pcaResult.p = p; % to save the structure p containing input parameters under the pcaResult
pcaResult.unitIdx = unitIdxPCA; % to save the index of units used in pca 

if p.Results.fullCrossVal 
    pcaResult.nonCrossValDataSaveDir = fullfile(p.Results.filePath, sprintf('mat_results/run%03d', pcaRunNumber)); % save the directory where the corresponding data are saved automatically by the neuralTraj.m
    pcaResult.crossValDataSaveDir = fullfile(p.Results.filePath, sprintf('mat_results/run%03d', pcaRunNumber+1)); % save the directory where the corresponding cross-validation data are saved automatically by the neuralTraj.m
else
    pcaResult.nonCrossValDataSaveDir = fullfile(p.Results.filePath, sprintf('mat_results/run%03d', pcaRunNumber)); % save the directory where the corresponding data are saved automatically by the neuralTraj.m
end

saveName = sprintf('%s_%dD_%02dmsBin', p.Results.saveNameTag, p.Results.pcaDim, p.Results.pcaBinSize );
save(saveName, 'pcaResult', 'pcaDat')

% Generate elementary plots for sanity check 
% plot3D(pcaResult.seqTrain, 'xorth', 'dimsToPlot', 1:3);
% plot(pcaResult.currentParams.C(:,1))
% plot(pcaResult.seqTrain(1).xsm(3,:))
% plot(pcaResult.seqTrain(1).xorth(5,:))
% plotEachDimVsTime(pcaResult.seqTrain, 'xsm', pcaResult.binWidth);
% plotEachDimVsTime(pcaResult.seqTrain, 'xorth', pcaResult.binWidth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_runPCA( filePath, fileName, saveNameTag, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function 
        
        %parse input, and extract name-value pairs for the main function
        % 'runPCA'
                
        default_frHighPass  = 1;  % low FR cut (3Hz)  
        default_frLowPass   = 60; % high FR cut (60Hz) % low pass might not be necessary 
        default_timeRange   = [-0.5e3 0.5e3]; % in ms (e.g. -500 to 500 ms relative to the event at time 0)
        default_pcaBinSize = 50; % bin size to be used for PCA (25ms)
        default_pcaDim     = 5;  % default dimensionality for PCA
        default_pcaKernSD  = 50; % default SD for kernel smoothing in PCA
        default_fullCrossVal= false; % default boolean to run full cross validation or not 
        default_xcorThresholdPer = 0.2; % default cross-correlation spike cooccurrence tolerance (i.e., 20% of the total spike count of the fewer spiking unit of the pair)
        default_crossValFolds = 4;  % the number of cross-validation fold
        
        p = inputParser; % create parser object
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileName'); % file name to open
        addRequired(p,'saveNameTag'); % save name to save pcaResult
        
        addParameter(p,'frHighPass', default_frHighPass)        
        addParameter(p,'frLowPass', default_frLowPass) % FR high cut might not be necessary 
        addParameter(p,'timeRange', default_timeRange)
        addParameter(p,'pcaBinSize', default_pcaBinSize)
        addParameter(p,'pcaDim', default_pcaDim)
        addParameter(p,'pcaKernSD', default_pcaKernSD)
        addParameter(p,'fullCrossVal', default_fullCrossVal)
        addParameter(p,'xcorThresholdPer', default_xcorThresholdPer)
        addParameter(p,'crossValFolds', default_crossValFolds)
        
        parse(p, filePath, fileName, saveNameTag, vargs{:})
        
    end

end

