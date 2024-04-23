function [corrRiseTS, corrFallTS, pulseTrainIdx] = detecteventbythreshold_interactive_timeout(filePath, timeseries, sampFreq, detectTimeout, saveName, varargin)
%This function detects event markers, and output rise and fall time points relative to threshold.
% timeseries: raw data from which events will be detected
% sampFreq: sampling rate (Hz)
% detectTimeout: detection timeout (ms)
close all;

p = parse_input_detectevent_interactive_timeout( filePath, timeseries, sampFreq, detectTimeout, saveName, varargin ); % parse input
stdFactor = p.Results.stdFactor;
earlyCutoff = p.Results.earlyCutoff; 
lateCutoff = p.Results.lateCutoff; 

% Define your default values for the GUI (uicontrol)
defaultStdFactor = num2str(stdFactor); % as an example
defaultEarlyCutoff = num2str(earlyCutoff); % as an example
defaultLateCutoff = num2str(lateCutoff); % as an example

userSatisfied = false;

while ~userSatisfied

    %figure; hold on; plot(timeseries)
    if p.Results.filterArtifact
        timeseries = timeseries-min(timeseries);
        [b, a] = butter(4, p.Results.lowpassCutoff/p.Results.nyquist, 'low'); % 4th order Butterworth high-pass
        filtered_timeseries = filter(b, a, timeseries);
        timeseries = timeseries-filtered_timeseries;
    end

    meanTS = mean(timeseries); % mean timeseries
    stdTS  = std(timeseries);  % std timeseries

    if isempty(p.Results.manualthresTS)
        thresTS = meanTS + stdFactor*stdTS; % detection threshold
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
    valRiseTS = valRiseTS(valRiseTS>earlyCutoff & valRiseTS<lateCutoff);
    valFallTS = valFallTS(valFallTS>earlyCutoff & valFallTS<lateCutoff);

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
                    end
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

    % ensure that the corrFallTS entries are larger than corrRiseTS row-by-row
    if (length(corrRiseTS) ~= length(corrFallTS)) || sum(corrFallTS-corrRiseTS<0)>=1

        corrFallTS_srt = nan(size(corrRiseTS, 1), size(corrRiseTS, 2));
        for jj = 1:length(corrRiseTS)
            tempCorrFallTSIdx = find(corrRiseTS(jj) < corrFallTS, 1, 'first');
            if ~isempty(tempCorrFallTSIdx)
                corrFallTS_srt(jj) = corrFallTS(tempCorrFallTSIdx);
            end
        end
        corrFallTS = corrFallTS_srt;
        corrValI = ~isnan(corrFallTS) & ~isnan(corrRiseTS);

        corrFallTS = corrFallTS(corrValI);
        corrRiseTS = corrRiseTS(corrValI);

    end

    if p.Results.cutoffShort
        ts_dur = corrFallTS - corrRiseTS;
        corrRiseTS = corrRiseTS(ts_dur>p.Results.short*sampFreq);
        corrFallTS = corrFallTS(ts_dur>p.Results.short*sampFreq);
    end

    if p.Results.chunkPulses
        chunkRiseTS = corrRiseTS(diff([0, corrRiseTS])>p.Results.chunkInterval/1000*sampFreq); % this should detect pulses beginning each train (train beginner)
        [~,~,pulseTrainIdx] = histcounts(corrRiseTS,[chunkRiseTS, corrRiseTS(end)]);

        if p.Results.removeOffTrainPulses % spot the pulses outside of the trains (usually single pulses) and get rid of them
            uniqueTrains = unique(pulseTrainIdx); % unique pulse train IDs
            pulsesInTrains = cell2mat(arrayfun(@(a) sum(pulseTrainIdx==a), uniqueTrains, 'UniformOutput', false));
            offTrainPulseIds = uniqueTrains(pulsesInTrains<=5);

            valPulseTrainIdx = ~ismember(pulseTrainIdx, offTrainPulseIds);
            corrRiseTS = corrRiseTS(valPulseTrainIdx);
            corrFallTS = corrFallTS(valPulseTrainIdx & corrFallTS > corrRiseTS(1));
            pulseTrainIdx = pulseTrainIdx(valPulseTrainIdx);
        end
    else
        pulseTrainIdx = nan;
    end

    if p.Results.findEarlyOnset
        t1 = round(p.Results.earlyOnsetWindow * sampFreq);
        corrRiseTS = corrRiseTS(corrRiseTS>t1); 
        earlyWin = arrayfun(@(a) timeseries(a-t1:a-1), corrRiseTS, 'un', 0);
        earlyThres = cellfun(@(a) (max(a)-min(a))/2+min(a), earlyWin, 'un', 0);
        earlyThresPt = cellfun(@(a, b) find(a>=b, 1, 'first'), earlyWin, earlyThres, 'un', 0);
        adjust = cell2mat(cellfun(@(a) t1-a, earlyThresPt, 'un', 0));
        corrRiseTS_late = corrRiseTS;
        corrRiseTS = corrRiseTS-adjust;
    end

    if exist(fullfile(filePath, 'Figure'), 'dir') ~= 7
        mkdir(fullfile(filePath, 'Figure'))
    end

    if p.Results.plotRez
        f = figure; hold on;
        plot(timeseries);
        plot(corrFallTS, ones(1, length(corrFallTS)) .* thresTS, '*r');
        if p.Results.findEarlyOnset
            plot(corrRiseTS, cell2mat(earlyThres), 'og');
            plot(corrRiseTS_late, ones(1, length(corrRiseTS)) .* thresTS, 'ob');
        else
            plot(corrRiseTS, ones(1, length(corrRiseTS)) .* thresTS, 'ob');
        end
        title([saveName sprintf(': %d events detected with stdFactor = %.1f', length(corrRiseTS), stdFactor)])
        %xlim([0 round(length(timeseries) / 3)])
        hold off;

        screenSize = get(0, 'ScreenSize');
        screenWidth = screenSize(3);
        screenHeight = screenSize(4);

        uf = figure('Visible', 'on', 'Position', [screenWidth / 2 - 150, screenHeight / 2 - 150, 300, 250], ...
            'MenuBar', 'none', 'Name', 'User Input', 'NumberTitle', 'off', 'WindowStyle', 'normal');

        uicontrol('Style', 'text', 'String', 'Do you accept the pulse detection result?', ...
            'Position', [15, 190, 270, 30], 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 12);

        % "Yes" button callback using a function handle
        yesCallback = @(src,event) dealWithYesButton(src);
        uicontrol('Style', 'pushbutton', 'String', 'Yes', 'Position', [45, 155, 100, 40], 'Callback', yesCallback); % Yes Button

        % "No" button callback
        stdFactorLabel = uicontrol('Style', 'text', 'String', 'Input the new StdFactor:', ...
            'Position', [15, 135, 270, 20], 'HorizontalAlignment', 'center', 'FontSize', 12, 'Visible', 'off');
        stdFactorField = uicontrol('Style', 'edit', 'String', defaultStdFactor, 'Position', [100, 110, 100, 25], 'Visible', 'off');
        
        timeCutoffLabel = uicontrol('Style', 'text', 'String', 'Input the early and late cutoffs:', ...
           'Position', [15, 80, 270, 20], 'HorizontalAlignment', 'center', 'FontSize', 12, 'Visible', 'off');
        earlyCutoffField = uicontrol('Style', 'edit', 'String', defaultEarlyCutoff, 'Position', [50, 55, 80, 25], 'Visible', 'off'); 
        lateCutoffField = uicontrol('Style', 'edit', 'String', defaultLateCutoff, 'Position', [180, 55, 80, 25], 'Visible', 'off');

        submitButton = uicontrol('Style', 'pushbutton', 'String', 'Submit', 'Position', [110, 30, 80, 25], ... % Submit Button
            'Callback', @(src,event) dealWithSubmitButton(src, stdFactorField, earlyCutoffField, lateCutoffField), 'Visible', 'off');

        noCallback = @(src,event) dealWithNoButton(src, stdFactorLabel, stdFactorField, submitButton, timeCutoffLabel, earlyCutoffField, lateCutoffField); 
        uicontrol('Style', 'pushbutton', 'String', 'No', 'Position', [155, 155, 100, 40], 'Callback', noCallback); % NO Button

        set(uf, 'UserData', struct('buttonPressed', 'none', 'stdFactorValue', NaN));
        uiwait(uf);
        userData = get(uf, 'UserData');

        if strcmp(userData.buttonPressed, 'yes')
            userSatisfied = true;
            if ishandle(uf)
                close(uf);
            end
        else
            stdFactor = userData.stdFactorValue;
            earlyCutoff = userData.earlyCutoff; 
            lateCutoff = userData.lateCutoff; 

            if ishandle(uf)
                close(uf);
            end
            if ishandle(f)
                close(f);
            end
            fprintf('Redetecting events with the new stdFactor (threshold)!\n')
        end
    end
end % end of while loop

%Save the figure to the 'Figure' directory in the specified filePath
%savefig(f, fullfile(filePath, 'Figure', [p.Results.saveName '.fig']));
print(f, fullfile(filePath, 'Figure', [p.Results.saveName '.dpdf']), '-dpdf', '-painters');
close(f)

%% nested helper function
    function p = parse_input_detectevent_interactive_timeout(filePath, timeseries, sampFreq, detectTimeout, saveName, vargs)
        % parse input, and extract name-value pairs
        default_stdFactor = 1;          % std factor
        default_plotRez = false;        % boolean for plotting
        default_chunkPulses = true;     % boolean for pulse chunking to get trial-by-trial pulses
        default_chunkInterval = 1000;   % interval by which chunking pulses (in ms)
        default_earlyCutoff = 1;        % Note that it's data index (points) not time in seconds
        default_lateCutoff = length(timeseries); % Note that it's data index (points) not time in seconds
        default_correctLongPulse = false; % correct for the possible long pulses, especially in the trEnd
        default_manualthresTS = [];     % manual threshold
        default_cutoffShort = false;
        default_short = 1;              % 1 sec
        default_findEarlyOnset = false; % try to detect the early onset point
        default_earlyOnsetWindow = 0.5; % detect the early onset point within -0.5 to 0 window with 0 being the detected point
        default_filterArtifact = false; % use low-pass filtered signal to remove the photodiode artifact
        default_lowpassCutoff = 5;      % 5Hz
        default_nyquist = 100;          % the frequency of the main signal to be detected (e.g., faceCam pulses are 200Hz)
        default_removeOffTrainPulses = false; % to remove lone artifact pulses outside of the train pulses

        p = inputParser; % create parser object
        addRequired(p,'filePath')
        addRequired(p,'timeseries')
        addRequired(p,'sampFreq')
        addRequired(p,'detectTimeout')
        addRequired(p,'saveName')
        addParameter(p,'stdFactor', default_stdFactor)
        addParameter(p,'plotRez', default_plotRez)
        addParameter(p,'chunkPulses', default_chunkPulses)
        addParameter(p,'chunkInterval', default_chunkInterval)
        addParameter(p,'earlyCutoff', default_earlyCutoff) % Note that it's data index (points) not time in seconds
        addParameter(p,'lateCutoff', default_lateCutoff) % Note that it's data index (points) not time in seconds
        addParameter(p,'correctLongPulse', default_correctLongPulse)
        addParameter(p,'manualthresTS', default_manualthresTS)
        addParameter(p,'cutoffShort', default_cutoffShort)
        addParameter(p,'short', default_short)
        addParameter(p,'findEarlyOnset', default_findEarlyOnset)
        addParameter(p,'earlyOnsetWindow', default_earlyOnsetWindow)
        addParameter(p,'filterArtifact', default_filterArtifact)
        addParameter(p,'lowpassCutoff', default_lowpassCutoff)
        addParameter(p,'nyquist', default_nyquist)
        addParameter(p,'removeOffTrainPulses', default_removeOffTrainPulses)

        parse(p, filePath, timeseries, sampFreq, detectTimeout, saveName, vargs{:})
    end

% Separate callback functions
    function dealWithYesButton(src)
        set(src.Parent, 'UserData', struct('buttonPressed', 'yes', 'stdFactorValue', NaN));
        uiresume(src.Parent);
    end

    function dealWithNoButton(src, stdFactorLabel, stdFactorField, submitButton, timeCutoffLabel, earlyCutoffField, lateCutoffField)
        set(stdFactorLabel, 'Visible', 'on');
        set(stdFactorField, 'Visible', 'on');
        set(submitButton, 'Visible', 'on');
        set(timeCutoffLabel, 'Visible', 'on');
        set(earlyCutoffField, 'Visible', 'on');
        set(lateCutoffField, 'Visible', 'on');
    end

    function dealWithSubmitButton(src, stdFactorField, earlyCutoffField, lateCutoffField)
        userData = get(src.Parent, 'UserData');
        userData.buttonPressed = 'submit';
        userData.stdFactorValue = str2double(get(stdFactorField, 'String'));
        userData.earlyCutoff = str2double(get(earlyCutoffField, 'String')); 
        userData.lateCutoff = str2double(get(lateCutoffField, 'String'));  
        set(src.Parent, 'UserData', userData);
        uiresume(src.Parent);
    end

   


end

