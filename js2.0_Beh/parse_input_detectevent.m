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