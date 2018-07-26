function [ binSpkCount ] = psthBinBehavNeural( neural, behav, binsize, psthWin, stableTrials )
%This function generates spike count vectors with the bin size and the size of the psth window specified by the user
% INPUT:
%   neural: neural spike times data as a cell (1-by-# of units)
%           Spike times are assumed to be in milliseconds (downsampled at 1kHz)
%   behav:  behavioral event time vector (it has to be a trial-by-1 vector) to which the spike times will be aligned
%   binsize: bin size (e.g. 20 ms)
%   psthWin: psth window centering around the event of interest e.g. [ 1e3, 1e3 ]
%   stableTrials: to or not to exclude the initial trials before stabilization of neural activity (default = -1, no exclusion of trials), or the # of initial trials for exclustion can be entered
%
% OUTPUT: binned spike count data as a structure (1-by-# of units)
%   binSpkCount: binSpkCount is a structure containing binned spike counts of
%   all units recorded from an animal aligned to the task event (variable) of interest
%   binSpkCount.SpkCountMat contains the binned spike count matrices  
%   binSpkCount.trAVG contains the mean spike count of each unit (across trials)
%   binSpkCount.trSEM contains the sem spike count of each unit (across trials)

% primary inputs
popData = neural;       % neural population data (cell) - spike times
alignBehVar = behav;    % the behavioral events to which neural data are aligned

% primary output
binSpkCount = struct;   % define binSpkCount as a structure

%% Set relevant parameters (bin edges, Gaussian kernel)
% specify bin edges
currParams.binSize  = binsize;    % store bin size
currParams.psthWin  = psthWin;    % store the psth window
currParams.binEdges = linspace(-psthWin(1,1), psthWin(1,1), sum(psthWin)/binsize+1); % binEdges

% set a gaussian kernel to be convolved with the spike train (delta function)
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 15;  % gaussian std
[currParams.filter.kernel]     = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

binSpkCount.currParams = currParams;  % store the relevant binning parameters

currEvt = alignBehVar(:,1); % current event vector

for u = 1:length(popData)   % # of units
    
    % get delta function for the whole spike train
    numAllSpks = length(popData{1,u});                       % # of spikes of the current unit
    delta = zeros(1,ceil(popData{1,u}(numAllSpks)));         % delta function of the spike train
    delta(1,ceil(popData{1,u})) = 1;                         % delta function of the spike train - put 1 at spike times
    tmpSmooth = conv(delta,currParams.filter.kernel,'same'); % convultion of the delta function with the Gaussian kernel
    
    tmpBinCountMat = [];       % temporary binned spike count matrix for the current unit
    
    for tr = 1:length(currEvt) % # of trial of the current event
        [ tmpBinCount ] = histcounts(popData{1,u}, currEvt(tr) + currParams.binEdges); % histogram bin counts
        tmpBinCountMat = [tmpBinCountMat; tmpBinCount]; % append the binned spike counts of the current trial to the current unit's matrix
    end
    
    %tmpBinCountMat = tmpBinCountMat.*(1000/binsize);     % to convert counts to Hz (optional)
    binSpkCount(u).SpkCountMat = sparse(tmpBinCountMat); % store the binned spike count matrix into the binSpkCount structure
    
    % get the mean and the std SC across trials
    if stableTrials==-1 % in this case, include all trials without excluding any trials in the beginning
        binSpkCount(u).trAVG        = mean(tmpBinCountMat,1); % mean spike counts averaged across trials
        binSpkCount(u).trSEM        = std(tmpBinCountMat,0,1) ./ sqrt(size(tmpBinCountMat,1)-1); % sem spike counts
    else                % in this case, exclude initial trials from mean and std calculation
        binSpkCount(u).trAVG        = mean(tmpBinCountMat(stableTrials:count,:),1); % mean spike counts averaged across trials
        binSpkCount(u).trSEM        = std(tmpBinCountMat(stableTrials:count,:),0,1) ./ sqrt(size(tmpBinCountMat(stableTrials:count,:),1)-1); % sem spike counts
    end
    
    %% normalization (Z-score): use the half of the left (prior to evt occurrence; time 0) window 
    avgSCz = mean(binSpkCount(u).trAVG(1,1:round(psthWin(1)/(binsize*2))), 2); % mean spike count of the bins in the left half of the psth window
    
    %  to bin the entire spike train (mean and std need to be compuated after binning)
    if alignBehVar(length(alignBehVar)) > numel(delta)
        currEntireBinEdges = linspace( alignBehVar(1), numel(delta), (numel(delta)-alignBehVar(1)+1)/binsize +1 );
    else
        currEntireBinEdges = linspace( alignBehVar(1), alignBehVar(length(alignBehVar)), (alignBehVar(length(alignBehVar))-alignBehVar(1)+1)/binsize +1 );
    end
    [ tmpEntireBinCount ] = histcounts(popData{1,u}, currEntireBinEdges); % binned spike count applied to the relevant part of the spike train delta function
    stdSCz = std(tmpEntireBinCount); % std across time bins to be used for the z-score normalization
    
    if stdSCz > 0.0001  % in case the std is large enough
        binSpkCount(u).SpkCountMatZ = (binSpkCount(u).trAVG-avgSCz)/stdSCz; % take the z score 
        binSpkCount(u).Zlogic       = 1; % z score logic is true
    else                % in case the std is too small
        disp(['STD=' num2str( stdSCz ) ' | cannot z-score']);   
        binSpkCount(u).SpkCountMatZ = binSpkCount(u).trAVG;     % do not take the z score
        binSpkCount(u).Zlogic       = 0; % z score logic is false
    end   
    
end % end of unit iteration

end % end of function

