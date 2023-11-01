function [figHandle, lickTimes, lickTrain] = lickRasterGrammLickTimes( psthWin, classLabel, manualX, varargin )
%lickRasterGrammLickTimes uses gramm package to generate a lick raster plot containing
% a single (or multiple) group(s) of lick trains.
%  Input: Variable number of lick train groups can be input as a cell containing
%   trial-by-trial lick trains.

%lengthargin = length(varargin); % varargin only counts the variable inputs
%not the required input (psthWin).

%binEdges1ms = linspace(-psthWin(1,1), psthWin(1,2), sum(psthWin)+1);  

class = [];
lickTimes = []; % lick time indices
lickTrain = []; % binned lick train

%cval={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I'}; % class labels (too many classes can be problematic!)
binSize = 50; % size of the bin to further bin the lick trains (i.e. 50 ms)
binX = linspace(0,sum(psthWin),round(sum(psthWin)/binSize)+1);      % bin for lick trains

% gaussian kernel to be convolved with the psths
gaussianSigma    = 1;  % gaussian std
[gaussianKernel] = TNC_CreateGaussian(gaussianSigma.*15,gaussianSigma,gaussianSigma.*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

for gr = 1:length(varargin) % increment groups
    
    % combine lickTimesCell across groups
    thisLickTimesCell = varargin{gr}; % lickTimesCell of the current group
    lickTimes = [lickTimes; thisLickTimesCell]; % append lickTimes
    
    % combine lickTrainCell across groups
    thisLickTrainCell = cell(length(thisLickTimesCell),1); % lickTrainCell of the current group
    for t = 1:length(thisLickTimesCell) % increment trials
        
        thisLickTrainCell{t,1} = histcounts(thisLickTimesCell{t,1},binX); % count licks across bins
        thisLickTrainCell{t,1} = conv(thisLickTrainCell{t,1}.*(1000/binSize),gaussianKernel,'same');  % conversion to Hz
        
        %thisLickTrainCell{t,1} = zeros(1,sum(psthWin));
        %thisLickTrainCell{t,1}(thisLickTimesCell{t,1}) = 1;
    end
    clearvars t
    
    lickTrain = [lickTrain; thisLickTrainCell]; % append lickTrains of different groups
    
    % put class label
    thisClass = zeros(length(varargin{gr}),1); % class label
    thisClass(:,1) = gr;         % label each class
    class = [class; thisClass];  % append class label
    
end
clearvars gr

c = classLabel(class); % class label
%c = cval(class); % class label

clear g

% lick rasters
g(1,1)=gramm('x',lickTimes,'color',c);
g(1,1).geom_raster();
%g(1,1).set_title('geom_raster()');
g(1,1).set_names('x','Time (ms)','y','Trials');

if ~isequal(psthWin, manualX)
    g(1,1).axe_property('xlim',[psthWin(1)-manualX(1) psthWin(1)-manualX(1)+sum(manualX)]); % manual xlim
end

% lick histogram (without errorbars)
% stat_bin includes binning, so a cell with discrete lick times can be input
% g(1,2)=gramm('x',lickTimes,'color',c);
% g(1,2).stat_bin('nbins',25,'geom','line');
% g(1,2).set_title('stat_bin()');

% mean +- sem (or other measures of variability)
xAxis=linspace(0,sum(psthWin),round(sum(psthWin)/binSize));
g(1,2)=gramm('x',xAxis,'y',lickTrain,'color',c); 
g(1,2).stat_summary('type','sem','setylim',true); % setylim true to scale the plot by the summarized data (not by the underlying data points) 
g(1,2).set_names('x','Time (ms)','y','FR (Hz)');

if ~isequal(psthWin, manualX)
    g(1,2).axe_property('xlim',[psthWin(1)-manualX(1) psthWin(1)-manualX(1)+sum(manualX)]); % manual xlim
    %g(1,2).set_title('stat_summary()');
end

figHandle = figure('Position',[100 100 800 350]);
g.draw();


end

