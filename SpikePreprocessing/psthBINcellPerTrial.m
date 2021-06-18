function [ unitTimeBOut, binSpkOut, binSpk, params ] = psthBINcellPerTrial( neural, behav, binsize, psthWin )
%This function generates spike count matrix for a single trial with the bin size and the size of the psth window specified by the user
% INPUT:
%   neural: neural spike times data as a cell (7-by-# of units)
%           Spike times are assumed to be in milliseconds (downsampled at 1kHz)
%   behav:   behavioral event time to which the spike times will be aligned (single trial is assumed)
%   binsize: bin size (e.g. 20 ms) - it must be noted that SpkCountMat only
%       is affected by the bin size. Other output components such as SpkGCountMat or SpkGCountMatBase are not affected by the binsize, they are currently computed on 1-ms bin data
%   psthWin: psth window centering around the event of interest e.g. [ 1e3, 1e3 ]
%
% OUTPUT: binned spike count data as a structure (1-by-# of units)
%   binSpk: binSpk is a structure containing binned spike counts of all units recorded from an animal aligned to the task event (variable) of interest
%   binSpk.SpkCountMat contains the binned spike count matrices

% primary output
binSpk    = struct;   % define binSpk as a structure (include all fields)
binSpkOut = struct;   % define binSpkOut as a structure (include seleceted fields only to output due to the memory issue)

%% Set relevant parameters (bin edges, Gaussian kernel)
% specify bin edges
params.binSize  = binsize;    % store bin size
params.psthWin  = psthWin;    % store the psth window
params.binEdges = linspace(-psthWin(1,1), psthWin(1,2), sum(psthWin)/binsize+1); % binEdges using the user-specified binsize
params.binEdges1ms = linspace(-psthWin(1,1), psthWin(1,2), sum(psthWin)+1);      % binEdges with 1ms binsize (to be used for detecting spike times)

currEvt = behav(1); % current event (single trial assumed)

for u = 1:size(neural,2)   % increment units
    % get delta function for the whole spike train
    delta = zeros(1,ceil(neural{1,u}(end))); % delta function of the whole spike train of the current unit
    delta(1,ceil(neural{1,u})) = 1; % delta function of the whole spike train - put 1 at spike times
    % get the binned spike count for the current evt (single trial assumed)
    if currEvt+psthWin(1,2) < length(delta) % in case the current trial's psth doesn't go out of the upper bound
        if currEvt-psthWin(1,1) > 0
            binSpk.SpkTimes{u,1}{1}     = find(histcounts(neural{1,u}, currEvt + params.binEdges1ms)==1);  % histogram bin counts with 1ms binsize
        else % in case the currEvt out of the lower bound: put NaNs
            binSpk.SpkTimes{u,1}{1}     = nan;  % histogram bin counts with 1ms binsize
        end
    else % in case the currEvt out of the upper bound: put NaNs
        binSpk.SpkTimes{u,1}{1}     = nan;  % histogram bin counts with 1ms binsize
    end
    clearvars tr
    % just to store information to the structure
    binSpk.fileInfo{u,1}    = strcat('unit',num2str(u)); % store the essential file information (animal id and date of data collection for each cell)
    binSpk.currEvt{u,1}     = currEvt;      % store the task events and the logical for valid trials (with spikes counted)
    binSpk.Site{u,1}        = neural{3,u};  % store the site, from which the current unit was recorded
    binSpk.isStr{u,1}       = neural{5,u};  % store the logic, indicating whether the unit's str or ctx
    
    if u == 1
        binSpk.params   = params;      % put params into the binSpk structure
    end
    
    binSpk.geometry{u,1} = neural{4,u}; % put geometry into the binSpk structure
    binSpk.meanWF{u,1}   = neural{6,u}; % put mean waveform into the binSpk structure
    
    %fprintf('processed unit %d\n', u) % report unit progression
    
    %% Get the unit-Time-Trial mat
    tmpUnitTrialTimeMat  = cell2mat(getSpkCntMatFromSpkTimes( binSpk.SpkTimes{u}, binSpk.params )); % get current unit's spikeCountMat (1ms bin)
    [binSpkOut.unitTimeB(u,:)] = bin1msSpkCountMat(tmpUnitTrialTimeMat,params.binSize,params.binSize);  % get current unit's binned spikeCountMat (e.g. 20ms bin) 
    
    %% Select which fields to output
    binSpkOut.SpkTimes{u,1}       = binSpk.SpkTimes{u,1};
    binSpkOut.fileInfo{u,1}       = binSpk.fileInfo{u,1};
    binSpkOut.currEvt{u,1}        = binSpk.currEvt{u,1};
    binSpkOut.Site{u,1}           = binSpk.Site{u,1};
    binSpkOut.geometry{u,1} = binSpk.geometry{u,1};
    binSpkOut.meanWF{u,1} = binSpk.meanWF{u,1};
    binSpkOut.isStr{u,1} = binSpk.isStr{u,1}; 
    %clearvars tmp*
end % end of unit iteration
    binSpkOut.params = binSpk.params;
    unitTimeBOut = sparse(binSpkOut.unitTimeB); % output as a separate variable
end % end of function

