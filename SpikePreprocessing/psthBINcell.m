function [ binSpkOut, binSpk, params ] = psthBINcell( fileInfo, regionInfo, neural, behav, baseEvt, binsize, psthWin, stableTrials, psthPlotFlag )
%This function generates spike count vectors with the bin size and the size of the psth window specified by the user
% INPUT:
%   neural: neural spike times data as a cell (5-by-# of units)
%           Spike times are assumed to be in milliseconds (downsampled at 1kHz)
%   behav:  behavioral event time vector (it has to be a trial-by-1 vector) to which the spike times will be aligned
%   baseEvt: behavioral event time vector (it has to be a trial-by-1 vector) to which the spike times will be aligned for a baseline period for z-score normalization
%   binsize: bin size (e.g. 20 ms) - it must be noted that SpkCountMat only
%       is affected by the bin size. Other output components such as SpkGCountMat or SpkGCountMatBase are not affected by the binsize, they are currently computed on 1-ms bin data
%   psthWin: psth window centering around the event of interest e.g. [ 1e3, 1e3 ]
%   stableTrials: to or not to exclude the initial trials before stabilization of neural activity (default = -1, no exclusion of trials), or the # of initial trials for exclustion can be entered
%
% OUTPUT: binned spike count data as a structure (1-by-# of units)
%   binSpk: binSpk is a structure containing binned spike counts of
%   all units recorded from an animal aligned to the task event (variable) of interest
%   binSpk.SpkCountMat contains the binned spike count matrices
%   binSpk.trAVG contains the mean spike count of each unit (across trials)
%   binSpk.trSEM contains the sem spike count of each unit (across trials)

% primary output
binSpk    = struct;   % define binSpk as a structure (include all fields)
binSpkOut = struct;   % define binSpkOut as a structure (include seleceted fields only to output due to the memory issue)

% just to figure out the dimension of the psth plot
row = ceil(sqrt(size(neural,2)))-1; % # of rows (# of features)
col = ceil(sqrt(size(neural,2)))+1; % # of columns (# of units)

%% Set relevant parameters (bin edges, Gaussian kernel)
% specify bin edges
params.binSize  = binsize;    % store bin size
params.psthWin  = psthWin;    % store the psth window
params.binEdges = linspace(-psthWin(1,1), psthWin(1,2), sum(psthWin)/binsize+1); % binEdges using the user-specified binsize
params.binEdges1ms = linspace(-psthWin(1,1), psthWin(1,2), sum(psthWin)+1);      % binEdges with 1ms binsize (to be used for detecting spike times)

% set a gaussian kernel to be convolved with the spike train (delta function)
params.smthParams.rise     = 1;   % gaussian kernel parameter
params.smthParams.decay    = 15;  % gaussian std
[params.filter.kernel]     = TNC_CreateGaussian(params.smthParams.decay.*15,params.smthParams.decay,params.smthParams.decay.*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

currEvt = behav(:,1); % current event vector

for u = 1:size(neural,2)   % increment units
    
    % get delta function for the whole spike train
    delta = zeros(1,ceil(neural{1,u}(end))); % delta function of the whole spike train of the current unit
    delta(1,ceil(neural{1,u})) = 1; % delta function of the whole spike train - put 1 at spike times
    tmpSmooth = conv(delta,params.filter.kernel,'same'); % convolution of the delta function with the Gaussian kernel
    validEvts = zeros(length(currEvt),1); % valid task events being used for spike counts
      
    valEvtCnt = 0; % valid event counter
    
    % get trial-by-trial binned spike counts
    for tr = 1:length(currEvt) % # of trial of the current event
            
        if currEvt(tr)+psthWin(1,2) < length(delta) % in case the current trial's psth doesn't go out of the upper bound
            if currEvt(tr)-psthWin(1,1) > 0
                valEvtCnt = valEvtCnt + 1; % valid event (event that occurred between the first and last spikes of each unit) count
                %binSpk.SpkCountMat{u,1}{tr,1}  = sparse(histcounts(neural{1,u}, currEvt(tr) + params.binEdges));      % histogram bin counts; 
                binSpk.SpkGCountMat{u,1}{tr,1} = tmpSmooth(1,currEvt(tr)-psthWin(1,1)+1:currEvt(tr)+psthWin(1,2));    % event bin counts from the Gaussian convolved delta function;        % histogram bin counts; 
                binSpk.SpkTimes{u,1}{tr,1}     = find(histcounts(neural{1,u}, currEvt(tr) + params.binEdges1ms)==1);  % histogram bin counts with 1ms binsize
                validEvts(tr,1) = true; % meaning that the spike counts were computed for the current trial               
            else % in case the currEvt out of the lower bound: put NaNs
                %binSpk.SpkCountMat{u,1}{tr,1}  = nan(1,length(params.binEdges)-1); % histogram bin counts; 
                binSpk.SpkGCountMat{u,1}{tr,1} = nan(1,length(params.binEdges)-1); % event bin counts from the Gaussian convolved delta function;        % histogram bin counts; 
                binSpk.SpkTimes{u,1}{tr,1}     = nan;  % histogram bin counts with 1ms binsize
            end
        else % in case the currEvt out of the upper bound: put NaNs 
            %binSpk.SpkCountMat{u,1}{tr,1}  = nan(1,length(params.binEdges)-1); % histogram bin counts; 
            binSpk.SpkGCountMat{u,1}{tr,1} = nan(1,length(params.binEdges)-1); % event bin counts from the Gaussian convolved delta function;        % histogram bin counts; 
            binSpk.SpkTimes{u,1}{tr,1}     = nan;  % histogram bin counts with 1ms binsize
        end
        
    end 
    clearvars tr
    
    % get trial-by-trial baseline binned spike counts
    tmpBaseBinGCountMat = zeros(length(baseEvt),round(sum(psthWin)/2));   % temporary binned spike count matrix (Gaussian convolved) for the baseline period
    for tr = 1:length(baseEvt)
        if baseEvt(tr) < length(delta)
            if baseEvt(tr,1)-round(sum(psthWin)/2) >= 0 % to make sure there's enough (left) room for baseline period sampling
                tmpBaseBinGCountMat(tr,:) = tmpSmooth(1,baseEvt(tr)-round(sum(psthWin)/2)+1:baseEvt(tr)); % baseline bin counts from the Gaussian convolved delta function
            else
                tmpBaseBinGCountMat(tr,:) = nan; % baseline bin counts from the Gaussian convolved delta function
            end
        end
    end
    clearvars tr
    
    % just to store information to the structure
    binSpk.fileInfo{u,1}    = strcat(fileInfo,regionInfo,'unit',num2str(u)); % store the essential file information (animal id and date of data collection for each cell)
    binSpk.SpkGCountMatBase{u,1} = tmpBaseBinGCountMat; % store the baseline binned spike count matrix into binSpk structure
    binSpk.currEvt{u,1}     = [currEvt,validEvts];      % store the task events and the logical for valid trials (with spikes counted)
    binSpk.Site{u,1}        = neural{3,u};              % store the site, from which the current unit was recorded
    binSpk.isStr{u,1}       = neural{5,u};              % store the logic, indicating whether the unit's str or ctx
    
    % get the mean and the std SC across trials (these stats are to be used for z-score calculation, thus use the gaussian convolved spike trains)
    tmpBinGCountMat = cell2mat(binSpk.SpkGCountMat{u,1});
    if stableTrials==-1 % in this case, include all trials without excluding any trials in the beginning      
        binSpk.trAVG{u,1}        = nanmean(tmpBinGCountMat,1); % mean spike counts
        binSpk.trSEM{u,1}        = nanstd(tmpBinGCountMat,0,1) ./ sqrt(size(tmpBinGCountMat,1)-1); % sem spike counts
        binSpk.trAVGbase{u,1}    = nanmean(tmpBaseBinGCountMat,1);   % mean spike counts averaged across trials for the baseline period
    else                % in this case, exclude initial trials from mean and std calculation
        binSpk.trAVG{u,1}        = nanmean(tmpBinGCountMat(stableTrials:end,:),1); % mean spike counts averaged across trials
        binSpk.trSEM{u,1}        = nanstd(tmpBinGCountMat(stableTrials:end,:),0,1) ./ sqrt(size(tmpBinGCountMat(stableTrials:end,:),1)-1); % sem spike counts
        binSpk.trAVGbase{u,1}    = nanmean(tmpBaseBinGCountMat(stableTrials:end,:),1);   % mean spike counts averaged across trials for the baseline period
    end
    
    %% normalization (Z-score): use the mean spike count of the baseline period
    avgSCzBase = nanmean(binSpk.trAVGbase{u,1}, 2);   % mean spike count of bins in the baseline period
    stdSCzBase = nanstd(binSpk.trAVGbase{u,1}, 0, 2); % std spike count of bins in the baseline period
    
    if stdSCzBase > 0.0001  % in case the std is large enough
        binSpk.SpkCountMatZ{u,1}  = (binSpk.trAVG{u,1}-avgSCzBase)./stdSCzBase; %stdSCz; % take the z score
        binSpk.SpkCountMatZe{u,1} = binSpk.trSEM{u,1}./stdSCzBase;
        binSpk.Zlogic{u,1}        = 1; % z score logic is true
        
        %% plot psth
        if psthPlotFlag % if the boolean for psth plot is true, draw psths
            subplot(row,col,u);
            hold on;
            unitname = num2str(u);
            plotName = strcat('u#', unitname); % subplot name
            plot([-psthWin(1) psthWin(2)],[0 0],'-','Color',[0.2 0.2 0.2])
            ylim([-2 10]); % ylim zscores
            
            if u == 1 % for the first unit label the axes
                ylabel('Norm. Firing Rate');
                xlabel('Time (ms)');
            end
            
            shadedErrorBar(-psthWin(1):psthWin(2)-1, binSpk.SpkCountMatZ{u,1}, binSpk.SpkCountMatZe{u,1},'k');  % get shadedErrorBar function
            
            line([0 0],[-2 10]);
            title(plotName);
        else    % do not draw psths
        end
        
    else % in case the std is too small
        disp(['STD=' num2str( stdSCzBase ) ' | cannot z-score']);
        binSpk.SpkCountMatZ{u,1}  = binSpk.trAVG{u,1};     % do not take the z score
        binSpk.SpkCountMatZe{u,1} = binSpk.trSEM{u,1};     % do not take the z score
        binSpk.Zlogic{u,1}        = 0; % z score logic is false
    end
    
    if u == 1
        binSpk.params   = params;      % put params into the binSpk structure
    end
    
    binSpk.geometry{u,1} = neural{4,u}; % put geometry into the binSpk structure
    binSpk.meanWF{u,1}   = neural{6,u}; % put mean waveform into the binSpk structure
    
    fprintf('processed unit %d\n', u) % report unit progression
    
    %% Get the unit-Time-Trial mat
    tmpUnitTrialTimeMat  = cell2mat(getSpkCntMatFromSpkTimes( binSpk.SpkTimes{u}, binSpk.params )); % get the current unit's spikeCountMat (trial-by-1msBin)
    binSpkOut.unitTimeTrial(u,:,:) = permute(tmpUnitTrialTimeMat,[3 2 1]); % permute to get the unitTimeTrial mat
    
    %% Select which fields to output
    %binSpkOut.SpkCountMat{u,1}    = binSpk.SpkCountMat{u,1};  % fileSize is too large 
    %binSpkOut.SpkGCountMat{u,1}   = binSpk.SpkGCountMat{u,1}; % fileSize is too large
    binSpkOut.SpkTimes{u,1}       = binSpk.SpkTimes{u,1};
    binSpkOut.fileInfo{u,1}       = binSpk.fileInfo{u,1};
    %binSpkOut.SpkGCountMatBase{u,1} = binSpk.SpkGCountMat{u,1};
    binSpkOut.currEvt{u,1}        = binSpk.currEvt{u,1};
    binSpkOut.Site{u,1}         = binSpk.Site{u,1};
    binSpkOut.trAVG{u,1}        = binSpk.trAVG{u,1};
    binSpkOut.trSEM{u,1}        = binSpk.trSEM{u,1};
    %binSpkOut.trAVGbase{u,1}    = binSpk.trAVGbase{u,1};
    binSpkOut.SpkCountMatZ{u,1} = binSpk.SpkCountMatZ{u,1};
    %binSpkOut.SpkCountMatZe{u,1} = binSpk.SpkCountMatZe{u,1};
    binSpkOut.Zlogic{u,1} = binSpk.Zlogic{u,1};
    binSpkOut.geometry{u,1} = binSpk.geometry{u,1};
    binSpkOut.meanWF{u,1} = binSpk.meanWF{u,1};
    binSpkOut.isStr{u,1} = binSpk.isStr{u,1}; 
    %clearvars tmp*
end % end of unit iteration
    binSpkOut.params = binSpk.params;
end % end of function

