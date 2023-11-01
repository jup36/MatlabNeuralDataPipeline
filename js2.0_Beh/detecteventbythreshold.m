function [corrRiseTS, corrFallTS, pulseTrainIdx] = detecteventbythreshold(timeseries, sampFreq, detectTimeout, varargin)
%This function detects event markers, and output rise and fall time points relative to threshold.
% timeseries: raw data from which events will be detected
% sampFreq: sampling rate (Hz)
% detectTimeout: detection timeout (ms)

p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, varargin ); % parse input
%p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, {} ); % use this line instead when running line-by-line
%p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, {'detectLater',trStartIdx(1),'correctLongPulse',true} ); % use this line instead when running line-by-line

if p.Results.filterArtifact
    [b, a] = butter(4, p.Results.lowpassCutoff/(p.Results.sigFreq/2), 'low'); % 4th order Butterworth high-pass
    timeseries = filter(b, a, timeseries);
end

meanTS = mean(abs(timeseries)); % mean timeseries
stdTS  = std(abs(timeseries));  % std timeseries

if isempty(p.Results.manualthresTS)
    thresTS = meanTS + p.Results.stdFactor*stdTS; % detection threshold
else
    thresTS = p.Results.manualthresTS; 
end

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

    % correct for long pulses that did not go low after a high
    addRiseTS = [];
    addFallTS = [];
    longPulseCnt = 0;
    if p.Results.correctLongPulse
        if length(valFallTS)==length(valRiseTS)
            pulseInterval = mode(valFallTS-valRiseTS(1:length(valFallTS))); % normal pulse width
            for ts = 1:length(valRiseTS)
                if ~isempty(find(valFallTS>valRiseTS(ts),1))
                    if valFallTS(find(valFallTS>valRiseTS(ts),1))-valRiseTS(ts)>pulseInterval*2 % this is a long pulse
                        longPulseCnt = longPulseCnt+1;
                        addFallTS(longPulseCnt) = valRiseTS(ts)+pulseInterval; % add a fall
                        addRiseTS(longPulseCnt) = valFallTS(find(valFallTS>valRiseTS(ts),1))-pulseInterval; % add a rise

                    else
                    end
                else
                end
            end
            corrRiseTS = sort([valRiseTS, addRiseTS]);
            corrFallTS = sort([valFallTS, addFallTS]);
        elseif length(valRiseTS)-length(valFallTS)==1 && sum(valFallTS-valRiseTS(1:end-1)>=0)==length(valFallTS)
            corrRiseTS= valRiseTS(1:end-1);
            corrFallTS= valFallTS;

        elseif length(valFallTS)~=length(valRiseTS)
            error('The # of pulse rises and falls do not match!')
        end
    else
        corrRiseTS= valRiseTS;
        corrFallTS= valFallTS;
    end

    % sanity check plot
    % plot(timeseries(1:25000*150)); hold on
    % plot(valRiseTS(valRiseTS<=25000*150),thresTS,'*B');
    % plot(valFallTS(valFallTS<=25000*150),thresTS,'*R'); hold off

    if p.Results.chunkPulses
        chunkRiseTS = corrRiseTS(diff([0, corrRiseTS])>p.Results.chunkInterval/1000*sampFreq); % this should detect pulses beginning each train (train beginner)
        [~,~,pulseTrainIdx] = histcounts(corrRiseTS,[chunkRiseTS, corrRiseTS(end)]);
    else
        pulseTrainIdx = nan;
    end

    if p.Results.plotRez
        figure; hold on;
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
        default_correctLongPulse = false; % correct for the possible long pulses, especially in the trEnd
        default_manualthresTS = [];     % manual threshold
        default_cutoffShort = false; 
        default_short = 1;              % 1 sec 
        default_findEarlyOnset = false; % try to detect the early onset point
        default_earlyOnsetWindow = 0.5; % detect the early onset point within -0.5 to 0 window with 0 being the detected point  
        default_filterArtifact = false; % use low-pass filtered signal to remove the photodiode artifact
        default_lowpassCutoff = 5;      % 5Hz
        default_sigFreq = 200;          % the frequency of the main signal to be detected (e.g., faceCam pulses are 200Hz) 


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
        addParameter(p,'correctLongPulse', default_correctLongPulse)
        addParameter(p,'manualthresTS', default_manualthresTS)
        addParameter(p,'cutoffShort', default_cutoffShort)
        addParameter(p,'short', default_short)
        addParameter(p,'findEarlyOnset', default_findEarlyOnset)
        addParameter(p,'earlyOnsetWindow', default_earlyOnsetWindow)
        addParameter(p,'filterArtifact', default_filterArtifact)
        addParameter(p,'lowpassCutoff', default_lowpassCutoff)
        addParameter(p,'sigFreq', default_sigFreq)


        parse(p,timeseries, sampFreq, detectTimeout, vargs{:})
end
end

