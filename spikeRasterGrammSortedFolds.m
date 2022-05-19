function [figHandle, spikeTimes, spikeTrain] = spikeRasterGrammSortedFolds( psthWin, manualX, foldDatCell )
%spikeRasterGrammSortedFolds uses gramm package to generate a spike raster plot containing
% a single (or multiple) group(s) of spike trains provided within a cell 'foldDatCell' .
%  Input: Variable number of spike train groups/folds can be input as a cell containing
%   trial-by-trial spike trains, but they MUST be bundled as a cell 'foldDatCell'.

class = [];
spikeTimes = []; % spike time indices
spikeTrain = []; % binned spike train

classLabel = cell(1,length(foldDatCell));  
foldIdFmt = 'Fold#%d'; 
for c = 1:length(classLabel)
    classLabel{1,c} = sprintf(foldIdFmt,c); 
end


%cval={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I'}; % class labels (too many classes can be problematic!)
binSize = 50; % size of the bin to further bin the spike trains (i.e. 50 ms)
binX = linspace(0,sum(psthWin),round(sum(psthWin)/binSize)+1);      % bin for spike trains

% gaussian kernel to be convolved with the psths
gaussianSigma    = 1;  % gaussian std
[gaussianKernel] = TNC_CreateGaussian(gaussianSigma.*10,gaussianSigma,gaussianSigma.*20,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

for fd = 1:length(foldDatCell) % increment folds
    
    % combine spikeTimesCell across groups
    thisSpikeTimesCell = foldDatCell{fd}; % spikeTimesCell of the current group
    spikeTimes = [spikeTimes; thisSpikeTimesCell]; % append spikeTimes
    
    % combine spikeTrainCell across groups
    thisSpikeTrainCell = cell(length(thisSpikeTimesCell),1); % spikeTrainCell of the current group
    for t = 1:length(thisSpikeTimesCell) % increment trials
        
        thisSpikeTrain1ms = zeros(1,sum(psthWin));
        if sum(~isnan(thisSpikeTimesCell{t,1}))>0
            thisSpikeTrain1ms(thisSpikeTimesCell{t,1}) = 1;
        end
        thisSpikeTrainCell{t,1} = bin1msSpkCountMat( thisSpikeTrain1ms, 50, 50, 'align', 'center' );        
        %thisSpikeTrainCell{t,1} = histcounts(thisSpikeTimesCell{t,1},binX); % count spikes across bins
        thisSpikeTrainCell{t,1} = conv(thisSpikeTrainCell{t,1}.*(1000/binSize),gaussianKernel,'same');  % conversion to Hz
    end
    clearvars t
    
    spikeTrain = [spikeTrain; thisSpikeTrainCell]; % append spikeTrains of different groups
    
    % put class label
    thisClass = zeros(length(foldDatCell{fd}),1); % class label
    thisClass(:,1) = fd;         % label each class
    class = [class; thisClass];  % append class label
    
end
clearvars fd

c = classLabel(class); % class label
%c = cval(class); % class label

clear g

% spike rasters
spikeTimesZeroCentered = cellfun( @(x) x-psthWin(1), spikeTimes, 'UniformOutput',false ); % zero center the spike times
g(1,1)=gramm('x',spikeTimesZeroCentered,'color',c);
g(1,1).geom_raster();
%g(1,1).set_title('geom_raster()');
g(1,1).set_names('x','Time (ms)','y','Trials');

if ~isequal(psthWin, manualX)
    g(1,1).axe_property('xlim',[psthWin(1)-manualX(1) psthWin(1)-manualX(1)+sum(manualX)]-psthWin(1)); % manual xlim
end

% spike histogram (without errorbars)
% stat_bin includes binning, so a cell with discrete spike times can be input
% g(1,2)=gramm('x',spikeTimes,'color',c);
% g(1,2).stat_bin('nbins',25,'geom','line');
% g(1,2).set_title('stat_bin()');

% mean +- sem (or other measures of variability)
xAxis=round(linspace(0,sum(psthWin),round(sum(psthWin)/binSize))-psthWin(1));
xAxisC = cell(length(spikeTrain),1); 
xAxisC = cellfun(@(a) xAxis, xAxisC, 'Un',0);  
g(1,2)=gramm('x',xAxisC,'y',spikeTrain,'color',c);
g(1,2).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points)
g(1,2).set_names('x','Time (ms)','y','FR (Hz)');

if ~isequal(psthWin, manualX)
    g(1,2).axe_property('xlim',[psthWin(1)-manualX(1) psthWin(1)-manualX(1)+sum(manualX)]-psthWin(1)); % manual xlim
    %g(1,2).axe_property('ylim',[4 14]);
    %g(1,2).set_title('stat_summary()');
end

figHandle = figure('Position',[100 100 800 350]);
g.draw();


end

