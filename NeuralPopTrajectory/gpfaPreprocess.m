function [ gpfaDat ] = runGPFA( filePath, fileName, timeRange, varargin )
%gpfaPreprocess prepares neural population data in a format to run gpfa. 
% The input structure should provide trial x time spikecount matrix. 
% The output structure 'gpfaDat' must be in the dimension of Trials with 2 fields (trialId, spikes). 
% Each entry (trial) of the field 'spikes' must be Neuron-by-Timebins matrix. 

cd(filePath)

p = parse_input_gpfa( filePath, fileName, varargin ); % parse input
% p = parse_input_gpfa( filePath, 'binSpkCountSTRIT01_121317pcaPSTHreach200ms.mat', {}); % when running line-by-line

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
if length(unique(cellfun(@length, S(:).SpkCountMat))) == 1 
   numbTrial = unique(cellfun(@length,S(:).SpkCountMat));
else
   error('There are units with different number of trials!!!')
end

valTimeBins = find( S.params{1}.binEdges1ms>=timeRange(1) & S.params{1}.binEdges1ms<timeRange(2) ); % find bins within the time range
unitMeanFR  = zeros(size(S.SpkCountMat,1),1); % unit mean FR rates within the valid time range to exclude low FR (e.g. < 2Hz) units
unitTimeTrial = zeros(size(S.SpkCountMat,1), length(valTimeBins), numbTrial); % unit x time(1ms bin) x trial matrix 

for u = 1:size(S.SpkCountMat,1) % increment units
    tmpUnitTrialTimeMat = full(cell2mat(S.SpkCountMat{u})); % trial x time spike count mat
    unitTimeTrial(u,:,:) = permute(tmpUnitTrialTimeMat(:,valTimeBins),[3 2 1]); %
    unitMeanFR(u) = sum(sum(unitTimeTrial(u,:,:)))/(numbTrial*length(valTimeBins))*1000; % get the mean FR(Hz) of the unit
end
clearvars u 

%% get rid of one unit of the pairs with high cross-correlation 
unitTimeAcrossTrials = reshape(unitTimeTrial, [size(unitTimeTrial,1), size(unitTimeTrial,2)*size(unitTimeTrial,3)]); % reshape the unitTimeTrial matrix to get unitTimeAcrossTrials (concatenate across trials)
xcorUnitIdx = getUnitOfHighXcorr( unitTimeAcrossTrials, 0.3 );
    
for tr = 1:numbTrial % increment trials
    gpfaDat(tr).trialId = tr; % trial ID
    gpfaDat(tr).spikes = unitTimeTrial(unitMeanFR > 3 & unitMeanFR < 50 & xcorUnitIdx,:,tr); % neuron x timeBin mat for the trial, to include pca-sorted units only, pc.pcUnitSort(pc.pcUnitSort(:,3)==1,1)
end
clearvars tr

%% Extract neural trajectories
result = neuralTraj(8, gpfaDat, 'datFormat', 'spikes', 'method', 'gpfa', 'xDim', 3, 'binWidth', 20, 'kernSDList', 50);

kernSD = 50;

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);

% Build 'seqTrainPlot' to plot specific trials of interest 
trialsToPlot = [1:3];  
counter = 0;    % counter for the trials to be plotted
for i = 1:length(seqTrain) % # of trials
    if sum(arrayfun(@(a) a==i, trialsToPlot))==1
        counter = counter+1;
        seqTrainPlot(counter).T = seqTrain(i).T;    % # of time windows
        seqTrainPlot(counter).trialId = seqTrain(i).trialId;
        seqTrainPlot(counter).xorth = seqTrain(i).xorth;
        seqTrainPlot(counter).xsm = seqTrain(i).xsm;
    else

    end
end

% Plot neural trajectories in 3D space
plot3D(seqTrainPlot, 'xorth', 'dimsToPlot', 1:3);

% Plot each dimension of neural trajectories versus time
plotEachDimVsTime(seqTrain, 'xorth', 3);

% Plot neural trajectories in 3D space
redTrials = b123mat(71:80,3);  % trials of which trajectories will be plotted in red
outlierTrial = b123mat(71,3);   % outlier trial of which trajectory will be plotted in a different color

plot3dJP(seqTrain, 'xsm', 'dimsToPlot', 1:3, 'redTrials', [100:110], 'outlierTrial', [1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_gpfa( filePath, fileName, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function 
        
        %parse input, and extract name-value pairs for the main function
        % 'PSTH_rasters'
                
        default_FRcut   = 2;  % default firing rate cut (e.g. 2 Hz) to exclude units with low FR
        default_binSize = 50; % default binSize (e.g. 50 ms)
        
        p = inputParser; % create parser object
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileName'); % file name to open
        
        addParameter(p,'FRcut', default_FRcut)        
        addParameter(p,'binSize', default_binSize)

        parse(p, filePath, fileName, vargs{:})
        
    end

end

