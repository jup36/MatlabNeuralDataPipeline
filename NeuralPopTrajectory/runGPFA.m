function [ gpfaResult, gpfaDat ] = runGPFA( filePath, fileName, saveNameTag, varargin )
%runGPFA prepares neural population data in a format to run gpfa. 
% The input structure should provide trial x time spikecount matrix. 
% The output structure 'gpfaDat' must be in the dimension of Trials with 2 fields (trialId, spikes). 
% Each entry (trial) of the field 'spikes' must be Neuron-by-Timebins matrix. 
% Modified on 9/6/18 to extract SpkCountMat from SpkTimes

p = parse_input_runGPFA( filePath, fileName, saveNameTag, varargin );
% p = parse_input_runGPFA( filePath, fileName, saveNameTag, {} ); % use this when running line-by-line

cd(p.Results.filePath)

gpfaFileFolder = fullfile(p.Results.filePath,'mat_results','run*'); % file folder where the gpfa results are saved  
gpfaRunNumber  = length(dir(gpfaFileFolder))+1; % get the right gpfaRunNumber to run neuralTraj.m, and save the result in the 'mat_results" folder under the current filePath. 
% this is just to avoid assigning 10 to the gpfaRunNumber, which for some
% reason would just keep skipping w/o running (BUG in the 'neuralTraj.m' )
if gpfaRunNumber == 10
    gpfaRunNumber = gpfaRunNumber+1;
end

load(p.Results.fileName, 'S', 'pc'); % load the raw data structure, and results of the pca analysis using 'pcaPSTH.m' 

% Check the existence of files
if exist('S','var') && ~exist('pc','var') 
    error('The variable pc is missing!!')
elseif ~exist('S','var') && exist('pc','var') 
    error('The variable S is missing!!')
elseif ~exist('S','var') && ~exist('pc','var') 
    error('Both S & pc are missing!!')
end

% Sanity check for the number of trials for each unit - check all units' spkCountMat have the same number of trials
if length(unique(cellfun(@length, S(:).SpkTimes))) == 1 
   numbTrial = unique(cellfun(@length,S(:).SpkTimes));
else
   error('There are units with different number of trials!!!')
end

valTimeBins = find( S.params.binEdges1ms>=p.Results.timeRange(1) & S.params.binEdges1ms<p.Results.timeRange(2) ); % find bins within the time range
unitMeanFR  = zeros(size(S.SpkTimes,1),1); % unit mean FR rates within the valid time range to exclude low and high FR (e.g. < 2Hz & >50Hz) units
unitTimeTrial = zeros(size(S.SpkTimes,1), length(valTimeBins), numbTrial); % unit x time(1ms bin) x trial matrix 

% organize the unitTimeTrial mat
for u = 1:size(S.SpkTimes,1) % increment units
    tmpUnitTrialTimeMat  = cell2mat(getSpkCntMatFromSpkTimes( S.SpkTimes{u}, S.params )); % get the current unit's spikeCountMat (trial-by-1msBin)
    %tmpUnitTrialTimeMat  = full(cell2mat(S.SpkCountMat{u})); % trial x time spike count mat
    unitTimeTrial(u,:,:) = permute(tmpUnitTrialTimeMat(:,valTimeBins),[3 2 1]); % permute to get the unitTimeTrial mat
    unitMeanFR(u) = sum(sum(unitTimeTrial(u,:,:)))/(numbTrial*length(valTimeBins))*1000; % get the mean FR(Hz) of the unit
end
clearvars u 

% get rid of one unit of the pairs with high cross-correlation (remove duplicate/split units)
unitTimeAcrossTrials = reshape(unitTimeTrial, [size(unitTimeTrial,1), size(unitTimeTrial,2)*size(unitTimeTrial,3)]); % reshape the unitTimeTrial matrix to get unitTimeAcrossTrials (concatenate across trials)
xcorUnitIdx = getUnitOfHighXcorr( unitTimeAcrossTrials, p.Results.xcorThresholdPer );
    
unitIdxGPFA = unitMeanFR > p.Results.frHighPass & unitMeanFR < p.Results.frLowPass & xcorUnitIdx; % units to be used for gpfa

if isfield(S,'nonNaNtrialId') % in case the nanNaNtrialId has been already sorted by the previous function 'pcaPSTHtwoPopulations.m'
    for tr = 1:length(S.nonNaNtrialId) % increment trials
        gpfaDat(tr).trialId = S.nonNaNtrialId(tr,1); % trial ID
        gpfaDat(tr).spikes = unitTimeTrial(unitIdxGPFA,:,S.nonNaNtrialId(tr,1)); % neuron x timeBin mat for the trial, to include pca-sorted units only, pc.pcUnitSort(pc.pcUnitSort(:,3)==1,1)
    end
else
    for tr = 1:numbTrial % increment trials
        gpfaDat(tr).trialId = tr; % trial ID
        gpfaDat(tr).spikes = unitTimeTrial(unitIdxGPFA,:,tr); % neuron x timeBin mat for the trial, to include pca-sorted units only, pc.pcUnitSort(pc.pcUnitSort(:,3)==1,1)
    end
    
end
clearvars tr

%% Extract neural trajectories 
gpfaResult = neuralTraj(gpfaRunNumber, gpfaDat, 'datFormat', 'spikes', 'method', 'pca', 'xDim', p.Results.gpfaDim, 'binWidth', p.Results.gpfaBinSize, 'kernSDList', p.Results.gpfaKernSD);

% Orthonormalize neural trajectories
[ gpfaResult.estParams, gpfaResult.seqTrain ] = postprocess(gpfaResult, 'kernSD', p.Results.gpfaKernSD); % orthonormalize

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.

% NOTES:
% - These function calls are computationally demanding.  Cross-validation 
%   takes a long time because a separate model has to be fit for each 
%   state dimensionality and each cross-validation fold.

if p.Results.fullCrossVal % cross-validation is not run by default
    
    for crossValXdim = [3 5 8]
        %neuralTraj(gpfaRunNumber+1, gpfaDat, 'method',  'pca', 'xDim', crossValXdim, 'numFolds', p.Results.crossValFolds);
        neuralTraj(gpfaRunNumber+1, gpfaDat, 'method', 'gpfa', 'xDim', crossValXdim, 'numFolds', p.Results.crossValFolds);
    end
    
    % Plot prediction error versus state dimensionality or versus kernSD.
        % Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
%     crossValKernSD = 50; % select kernSD for two-stage methods
%     plotPredErrorVsDim(gpfaRunNumber+1, crossValKernSD);
%     crossValXDim = 5; % select state dimensionality
%     plotPredErrorVsKernSD(gpfaRunNumber+1, crossValXDim);
    
end

%% Save gpfaResult separately (in addition to the automatic save by 'neuralTraj')
gpfaResult.p = p; % to save the structure p containing input parameters under the gpfaResult
gpfaResult.unitIdx = unitIdxGPFA; % to save the index of units used in gpfa 

if p.Results.fullCrossVal 
    gpfaResult.nonCrossValDataSaveDir = fullfile(p.Results.filePath, sprintf('mat_results/run%03d', gpfaRunNumber)); % save the directory where the corresponding data are saved automatically by the neuralTraj.m
    gpfaResult.crossValDataSaveDir = fullfile(p.Results.filePath, sprintf('mat_results/run%03d', gpfaRunNumber+1)); % save the directory where the corresponding cross-validation data are saved automatically by the neuralTraj.m
else
    gpfaResult.nonCrossValDataSaveDir = fullfile(p.Results.filePath, sprintf('mat_results/run%03d', gpfaRunNumber)); % save the directory where the corresponding data are saved automatically by the neuralTraj.m
end

saveName = sprintf('%s_%dD_%02dmsBin', p.Results.saveNameTag, p.Results.gpfaDim, p.Results.gpfaBinSize );
save(saveName, 'gpfaResult', 'gpfaDat')

% Generate elementary plots for sanity check 
% plot3D(gpfaResult.seqTrain, 'xorth', 'dimsToPlot', 1:3);
% plot(gpfaResult.currentParams.C(:,1))
% plot(gpfaResult.seqTrain(1).xsm(3,:))
% plot(gpfaResult.seqTrain(1).xorth(5,:))
% plotEachDimVsTime(gpfaResult.seqTrain, 'xsm', gpfaResult.binWidth);
% plotEachDimVsTime(gpfaResult.seqTrain, 'xorth', gpfaResult.binWidth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_runGPFA( filePath, fileName, saveNameTag, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function 
        
        %parse input, and extract name-value pairs for the main function
        % 'runGPFA'
                
        default_frHighPass  = 3;  % low FR cut (3Hz)  
        default_frLowPass   = 60; % high FR cut (60Hz) % low pass might not be necessary 
        default_timeRange   = [-0.5e3 0.5e3]; % in ms (e.g. -500 to 500 ms relative to the event at time 0)
        default_gpfaBinSize = 50; % bin size to be used for GPFA (25ms)
        default_gpfaDim     = 5;  % default dimensionality for GPFA
        default_gpfaKernSD  = 50; % default SD for kernel smoothing in GPFA
        default_fullCrossVal= false; % default boolean to run full cross validation or not 
        default_xcorThresholdPer = 0.2; % default cross-correlation spike cooccurrence tolerance (i.e., 20% of the total spike count of the fewer spiking unit of the pair)
        default_crossValFolds = 4;  % the number of cross-validation fold
        
        p = inputParser; % create parser object
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileName'); % file name to open
        addRequired(p,'saveNameTag'); % save name to save gpfaResult
        
        addParameter(p,'frHighPass', default_frHighPass)        
        addParameter(p,'frLowPass', default_frLowPass) % FR high cut might not be necessary 
        addParameter(p,'timeRange', default_timeRange)
        addParameter(p,'gpfaBinSize', default_gpfaBinSize)
        addParameter(p,'gpfaDim', default_gpfaDim)
        addParameter(p,'gpfaKernSD', default_gpfaKernSD)
        addParameter(p,'fullCrossVal', default_fullCrossVal)
        addParameter(p,'xcorThresholdPer', default_xcorThresholdPer)
        addParameter(p,'crossValFolds', default_crossValFolds)
        
        parse(p, filePath, fileName, saveNameTag, vargs{:})
        
    end

end

