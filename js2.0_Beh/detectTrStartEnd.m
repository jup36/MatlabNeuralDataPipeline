function [trStartIdx, trEndIdx] = detectTrStartEnd(trStartTimeseries, trEndTimeseries, sampFreq, detectTimeout, varargin)
%This function detects event markers, and output rise and fall time points relative to threshold.
% timeseries: raw data from which events will be detected
% sampFreq: sampling rate (Hz)
% detectTimeout: detection timeout (ms)

p = parse_input_detectTrStartEnd( trStartTimeseries, trEndTimeseries, sampFreq, detectTimeout, varargin ); % parse input
%p = parse_input_detectTrStartEnd( trStartTimeseries, trEndTimeseries, sampFreq, detectTimeout, {} ); % use this line instead when running line-by-line

[trStartRiseTS, trStartFallTS] = detecteventbythreshold(trStartTimeseries, sampFreq, detectTimeout, varargin); 
[trEndRiseTS, trEndFallTS] = detecteventbythreshold(); 




if p.Results.plotRez
    hold on; 
    plot(timeseries); 
    plot(valRiseTS,thresTS,'ob'); 
    plot(valFallTS,thresTS,'*r'); hold off
end

    %% nested helper function
    function p = parse_input_detectTrStartEnd( trStartTimeseries, trEndTimeseries, sampFreq, detectTimeout, vargs )
        % parse input, and extract name-value pairs
        default_stdFactor = 1;          % std factor
        default_plotRez = false;        % boolean for plotting
        default_chunkPulses = true;     % boolean for pulse chunking to get trial-by-trial pulses
        default_chunkInterval = 1000;   % interval by which chunking pulses (in ms)
        default_detectLater = 1;        % detect events later than a certain timepoint to prevent detection of premature events
        default_detectEarlier = length(timeseries); % detect events earlier than a certain timepoint to preclude late events
        
        p = inputParser; % create parser object
        addRequired(p,'trStartTimeseries')
        addRequired(p,'trEndTimeseries')
        addRequired(p,'sampFreq')
        addRequired(p,'detectTimeout')
        addParameter(p,'stdFactor', default_stdFactor)
        addParameter(p,'plotRez', default_plotRez)
        addParameter(p,'chunkPulses', default_chunkPulses)
        addParameter(p,'chunkInterval', default_chunkInterval)
        addParameter(p,'detectLater', default_detectLater)
        addParameter(p,'detectEarlier', default_detectEarlier)
        %addParameter(p,'detectFall', default_detectFall)
        
        parse(p,trStartTimeseries, trEndTimeseries, sampFreq, detectTimeout, vargs{:})
    end

end

