function [corrRiseTS, corrFallTS, pulseTrainIdx] = detecteventbythreshold_noAbs(timeseries, sampFreq, detectTimeout, varargin)
%This function detects event markers, and output rise and fall time points relative to threshold.
% timeseries: raw data from which events will be detected
% sampFreq: sampling rate (Hz)
% detectTimeout: detection timeout (ms)

p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, varargin ); % parse input
%p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, {} ); % use this line instead when running line-by-line
%p = parse_input_detectevent( timeseries, sampFreq, detectTimeout, {'detectLater',trStartIdx(1),'correctLongPulse',true} ); % use this line instead when running line-by-line

meanTS = mean(timeseries); % mean timeseries
stdTS  = std(timeseries);  % std timeseries

if isempty(p.Results.manualthresTS)
    thresTS = meanTS + p.Results.stdFactor*stdTS; % detection threshold
else
    thresTS = p.Results.manualthresTS;
end

detectInterval = sampFreq/1000*detectTimeout; % in ms

riseTS  = find(timeseries>thresTS);
valRiseTS = riseTS(diff([0, riseTS])>detectInterval);

% detect fall time points - just detect the threshold crossing after
% flipping the array
if isrow(timeseries)
    flipFallTS = find(fliplr(timeseries)>thresTS);
    flipTS = fliplr(1:length(timeseries));
else
    flipFallTS = find(flipud(timeseries)>thresTS);
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

if p.Results.cutoffShort
   ts_dur = corrFallTS - corrRiseTS; 
   corrRiseTS = corrRiseTS(ts_dur>p.Results.short*sampFreq); 
   corrFallTS = corrFallTS(ts_dur>p.Results.short*sampFreq); 
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

if p.Results.findEarlyOnset 
    t1 = round(p.Results.earlyOnsetWindow * sampFreq); 
    earlyWin = arrayfun(@(a) timeseries(a-t1:a-1), corrRiseTS, 'un', 0); 
    earlyThres = cellfun(@(a) (max(a)-min(a))/2+min(a), earlyWin, 'un', 0); 
    earlyThresPt = cellfun(@(a, b) find(a>=b, 1, 'first'), earlyWin, earlyThres, 'un', 0); 
    adjust = cell2mat(cellfun(@(a) t1-a, earlyThresPt, 'un', 0)); 
    corrRiseTS_late = corrRiseTS; 
    corrRiseTS = corrRiseTS-adjust;
end

if p.Results.plotRez
    figure; hold on;
    plot(timeseries);
    plot(corrFallTS,ones(1, length(corrFallTS)).*thresTS,'*r'); 
    if p.Results.findEarlyOnset 
        plot(corrRiseTS,cell2mat(earlyThres),'og');
        plot(corrRiseTS_late,ones(1, length(corrRiseTS)).*thresTS,'ob');
    else
        plot(corrRiseTS,ones(1, length(corrRiseTS)).*thresTS,'ob');
    end 
    hold off
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

        parse(p,timeseries, sampFreq, detectTimeout, vargs{:})
    end
end

