function [valRiseTS, valFallTS, pulseTrainIdx] = detecteventbythreshold(timeseries, sampFreq, detectTimeout, varargin)
%This function detects event markers, and output rise and fall time points relative to threshold.
% timeseries: raw data from which events will be detected
% sampFreq: sampling rate (Hz)
% detectTimeout: detection timeout (ms)

p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, varargin ); % parse input
%p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, {} ); % use this line instead when running line-by-line

meanTS = mean(abs(timeseries)); % mean timeseries
stdTS  = std(abs(timeseries));  % std timeseries
thresTS = meanTS + p.Results.stdFactor*stdTS; % detection threshold
detectInterval = sampFreq/1000*detectTimeout; % in ms

riseTS  = find(abs(timeseries)>thresTS);
valRiseTS = riseTS(diff([0, riseTS])>detectInterval);

% detect fall time points - just detect the threshold crossing after
% flipping the array
if isrow(timeseries)
    flipFallTS = find(fliplr(abs(timeseries))>thresTS);
    flipTS = fliplr(1:length(timeseries));
else
    flipFallTS = find(flipud(abs(timeseries))>thresTS);
    flipTS = flipud(1:length(timeseries));
end
valFallTS = sort(flipTS(flipFallTS(diff([0,flipFallTS])>detectInterval)),'ascend');

% check to exclude any early or late timepoints
valRiseTS = valRiseTS(valRiseTS>p.Results.detectLater & valRiseTS<p.Results.detectEarlier); 
valFallTS = valFallTS(valFallTS>p.Results.detectLater & valFallTS<p.Results.detectEarlier); 

% sanity check plot
% plot(timeseries(1:25000*150)); hold on
% plot(valRiseTS(valRiseTS<=25000*150),thresTS,'*B');
% plot(valFallTS(valFallTS<=25000*150),thresTS,'*R'); hold off

if p.Results.chunkPulses
    chunkRiseTS = riseTS(diff([0, riseTS])>p.Results.chunkInterval/1000*sampFreq); % this should detect pulses beginning each train (train beginner)
    [~,pulseTrainIdx] = min(abs(repmat(valRiseTS,length(chunkRiseTS),1)-chunkRiseTS'),[],1); % get the chunk id for each pulse by finding the closest train beginner
else
    pulseTrainIdx = nan;
end

if p.Results.plotRez
    hold on; 
    plot(timeseries); 
    plot(valRiseTS,thresTS,'ob'); 
    plot(valFallTS,thresTS,'*r'); hold off
end

    %% nested helper function
    function p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, vargs )
        % parse input, and extract name-value pairs
        default_stdFactor = 1;          % std factor
        default_plotRez = false;        % boolean for plotting
        default_chunkPulses = true;     % boolean for pulse chunking to get trial-by-trial pulses
        default_chunkInterval = 1000;   % interval by which chunking pulses (in ms)
        default_detectLater = 1;        % detect events later than a certain timepoint to prevent detection of premature events
        default_detectEarlier = length(timeseries); % detect events earlier than a certain timepoint to preclude late events
        
        p = inputParser; % create parser object
        addRequired(p,'timeseries')
        addRequired(p,'sampFreq')
        addRequired(p,'detectTimeout')
        addParameter(p,'stdFactor', default_stdFactor)
        addParameter(p,'plotRez', default_plotRez)
        addParameter(p,'chunkPulses', default_chunkPulses)
        addParameter(p,'chunkInterval', default_chunkInterval)
        addParameter(p,'detectLater', default_detectLater)
        addParameter(p,'detectEarlier', default_detectEarlier)
        %addParameter(p,'detectFall', default_detectFall)
        
        parse(p,timeseries, sampFreq, detectTimeout, vargs{:})
    end

end

