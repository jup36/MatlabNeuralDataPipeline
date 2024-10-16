function parse_tbyt_auditory_gng_behavior_auto(filePath, varargin)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_041624';
% varargin = {'preToneWin', 1, 'postToneWin', 4}

p = parseInput_tbyt_auditory_gng_behavior(filePath, varargin);
% p = parseInput_tbyt_auditory_gng_behavior(filePath, {'preToneWin', 1, 'postToneWin', 4}); 

%% locate the nidq folder
filePath_nidq = GrabFiles_sort_trials('_g', 0, {filePath});
if isempty(filePath_nidq)
    filePath_nidq = uigetdir(filePath, 'Select the nidq folder');
end

if length(filePath_nidq) > 1   
    warning("More than 1 nidq folders found! Will run with the first one!")
end
filePath_nidq = filePath_nidq{1}; % just take the path
 
%% Load evtInS or create one
evtInS = timestamp_behav_events_auto(filePath_nidq, false, 'cmosTrig', 'faceCam', 'speaker', ...
    'treadMill', 'topRed', 'sideGreen', 'water', 'airpuff', ...
    'blueLED', 'limeLED', 'greenLED', 'redLED', 'lick', 'manualWater'); % behavioral events

%% Load stimopts
filePath_stim = GrabFiles_sort_trials('_stimInfo', 0, {filePath});
if isempty(filePath_stim)
    [stimInfoName, stimInfoDir] = uigetfile(filePath, 'Select the stimInfo file!');
    load(fullfile(stimInfoDir, stimInfoName), 'stimopts')
end

if length(filePath_stim) > 1
    warning("More than 1 stimInfo files found! Will run with the first one!")
end

load(fullfile(filePath_stim{1}), 'stimopts'); 

%% Identify and sort events
%assert(size(evtInS.speaker, 1)==length(stimopts.stim_type), 'The number of recorded events and the number of events executed do not match!')
% Get the sizes
numRecordedEvents = size(evtInS.speaker, 1);
numExecutedEvents = length(stimopts.stim_type);

if numRecordedEvents ~= numExecutedEvents
    if numRecordedEvents > numExecutedEvents
        error('The number of recorded events exceeds the number of events executed!');
    else
        fprintf('The number of recorded events is less than the number of events executed.\n');
        userDecision = input('Do you want to proceed by cropping the executed events to match the recorded events? Y/N: ', 's');
        if strcmpi(userDecision, 'Y')
            % Crop stimopts.stim_type to match the length of evtInS.speaker
            stimopts.stim_type = stimopts.stim_type(1:numRecordedEvents);
            stimopts.rewarded_stim = stimopts.rewarded_stim(1:numRecordedEvents); 
            stimopts.punished_stim = stimopts.punished_stim(1:numRecordedEvents); 
            stimopts.outcome_positive_stim = stimopts.outcome_positive_stim(1:sum(stimopts.stim_type==stimopts.positive_stim)); 
            stimopts.outcome_negative_stim = stimopts.outcome_negative_stim(1:sum(stimopts.stim_type==stimopts.negative_stim)); 

            fprintf('Proceeding with cropped stimopts!\n');
        else
            error('User opted not to proceed. Exiting.');
        end
    end
else
    fprintf('The number of recorded and executed events match.\n');
end

evts = evtInS.speaker;
evts(:, 3) = stimopts.stim_type;

%% tbytDat
numbTr = size(evts, 1);
tbytDat = struct;
% event On
evtOnC = num2cell(evts(:, 1));
[tbytDat(1:numbTr).evtOn] = deal(evtOnC{:});
% event Off
evtOffC = num2cell(evts(:, 2));
[tbytDat(1:numbTr).evtOff] = deal(evtOffC{:});
% event type
evtTypeC = num2cell(evts(:, 3));
[tbytDat(1:numbTr).evtType] = deal(evtTypeC{:}); % 1: common auditory, 2: uncommon auditory

% define functions to detect relevant behavioral and task events
%detectTones = @(a, b) a <= evtInS.speaker(:, 1) & b >= evtInS.speaker(:, 1); % to detect tones
detectFaceCamPulses = @(a, b) evtInS.faceCam(a <= evtInS.faceCam(:, 1) & b >= evtInS.faceCam(:, 1)); % to detect faceCam exposure pulses
detectCmosGreenPulses = @(a, b) a <= evtInS.sideGreen(:, 1) & b >= evtInS.sideGreen(:, 1); % to detect cmos exposure pulses
detectCmosRedPulses = @(a, b) a <= evtInS.topRed(:, 1) & b >= evtInS.topRed(:, 1); % to detect cmos exposure pulses
detectBlueLEDpulses = @(a, b) a <= evtInS.blueLED(:, 1) & b >= evtInS.blueLED(:, 1); % to detect blue LED pulses
detectGreenLEDpulses = @(a, b) a <= evtInS.greenLED(:, 1) & b >= evtInS.greenLED(:, 1); % to detect green LED pulses
detectLimeLEDpulses = @(a, b) a <= evtInS.limeLED(:, 1) & b >= evtInS.limeLED(:, 1); % to detect lime LED pulses
detectRedLEDpulses = @(a, b) a <= evtInS.redLED(:, 1) & b >= evtInS.redLED(:, 1); % to detect red LED pulses
%detectLicks = @(a, b) evtInS.lick(a <= evtInS.lick & b >= evtInS.lick); % to detect peri-event licks
%detectWater = @(a, b) evtInS.water(a <= evtInS.water(:, 1) & b >= evtInS.water(:, 1)); % to detect water deliveries
%detectAirpuff = @(a, b) evtInS.airpuff(a <= evtInS.airpuff(:, 1) & b >= evtInS.airpuff(:, 1)); % to detect water deliveries

% assign block ID and CMOS exposure pulses
for tt = 1:length(tbytDat)
    % auditory stimuli trial
    % cmos pulses
    tempSideGreenI = detectCmosGreenPulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % cmos pulses
    tempTopRedI = detectCmosRedPulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % cmos pulses
    tempBlueLEDI = detectBlueLEDpulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % cmos pulses
    tempGreenLEDI = detectGreenLEDpulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % cmos pulses
    tempLimeLEDI = detectLimeLEDpulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % cmos pulses
    tempRedLEDI = detectRedLEDpulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % cmos pulses

    %% organize LED pulses
    % blue LED
    if sum(tempBlueLEDI)>0
        tbytDat(tt).blueLED = evtInS.blueLED(tempBlueLEDI, 1); % cmos pulses train Id
        tempBlueLEDTrainI = unique(evtInS.blueLED(tempBlueLEDI, 2));
        tbytDat(tt).blueLEDTrainI = tempBlueLEDTrainI;
        tempBlueLEDPulses = evtInS.blueLED(evtInS.blueLED(:, 2)==tempBlueLEDTrainI, 1);
        tempFirstPulseInBlueLEDTrain = find(ismember(tempBlueLEDPulses, tbytDat(tt).blueLED), 1, 'first');
        tbytDat(tt).blueLEDPulsesOfTrain = {tempFirstPulseInBlueLEDTrain, tempFirstPulseInBlueLEDTrain+numel(tbytDat(tt).blueLED)-1};
        if length(tbytDat(tt).blueLED) < (tbytDat(tt).evtOff-tbytDat(tt).evtOn + p.Results.preToneWin + p.Results.postToneWin - 1)*10 % 10Hz (40/4 LEDs)
            tbytDat(tt).blueLED = []; % omit insufficient number of frames
            tbytDat(tt).blueLEDTrainI = [];
            tbytDat(tt).blueLEDPulsesOfTrain = [];
        end
    end
    % green LED
    if sum(tempGreenLEDI)>0
        tbytDat(tt).greenLED = evtInS.greenLED(tempGreenLEDI, 1); % cmos pulses train Id
        tempGreenLEDTrainI = unique(evtInS.greenLED(tempGreenLEDI, 2));
        tbytDat(tt).greenLEDTrainI = tempGreenLEDTrainI;
        tempGreenLEDPulses = evtInS.greenLED(evtInS.greenLED(:, 2)==tempGreenLEDTrainI, 1);
        tempFirstPulseInGreenLEDTrain = find(ismember(tempGreenLEDPulses, tbytDat(tt).greenLED), 1, 'first');
        tbytDat(tt).greenLEDPulsesOfTrain = {tempFirstPulseInGreenLEDTrain, tempFirstPulseInGreenLEDTrain+numel(tbytDat(tt).greenLED)-1};
        if length(tbytDat(tt).greenLED) < (tbytDat(tt).evtOff-tbytDat(tt).evtOn + p.Results.preToneWin + p.Results.postToneWin - 1)*10 % 10Hz (40/4 LEDs)
            tbytDat(tt).greenLED = []; % omit insufficient number of frames
            tbytDat(tt).greenLEDTrainI = [];
            tbytDat(tt).greenLEDPulsesOfTrain = [];
        end
    end
    % lime LED
    if sum(tempLimeLEDI)>0
        tbytDat(tt).limeLED = evtInS.limeLED(tempLimeLEDI, 1); % cmos pulses train Id
        tempLimeLEDTrainI = unique(evtInS.limeLED(tempLimeLEDI, 2));
        tbytDat(tt).limeLEDTrainI = tempLimeLEDTrainI;
        tempLimeLEDPulses = evtInS.limeLED(evtInS.limeLED(:, 2)==tempLimeLEDTrainI, 1);
        tempFirstPulseInLimeLEDTrain = find(ismember(tempLimeLEDPulses, tbytDat(tt).limeLED), 1, 'first');
        tbytDat(tt).limeLEDPulsesOfTrain = {tempFirstPulseInLimeLEDTrain, tempFirstPulseInLimeLEDTrain+numel(tbytDat(tt).limeLED)-1};
        if length(tbytDat(tt).limeLED) < (tbytDat(tt).evtOff-tbytDat(tt).evtOn + p.Results.preToneWin + p.Results.postToneWin - 1)*10 % 10Hz (40/4 LEDs)
            tbytDat(tt).limeLED = []; % omit insufficient number of frames
            tbytDat(tt).limeLEDTrainI = [];
            tbytDat(tt).limeLEDPulsesOfTrain = [];
        end
    end
    % red LED
    if sum(tempRedLEDI)>0
        tbytDat(tt).redLED = evtInS.redLED(tempRedLEDI, 1); % cmos pulses train Id
        tempRedLEDTrainI = unique(evtInS.redLED(tempRedLEDI, 2));
        tbytDat(tt).redLEDTrainI = tempRedLEDTrainI;
        tempRedLEDPulses = evtInS.redLED(evtInS.redLED(:, 2)==tempRedLEDTrainI, 1);
        tempFirstPulseInRedLEDTrain = find(ismember(tempRedLEDPulses, tbytDat(tt).redLED), 1, 'first');
        tbytDat(tt).redLEDPulsesOfTrain = {tempFirstPulseInRedLEDTrain, tempFirstPulseInRedLEDTrain+numel(tbytDat(tt).redLED)-1};
        if length(tbytDat(tt).redLED) < (tbytDat(tt).evtOff-tbytDat(tt).evtOn + p.Results.preToneWin + p.Results.postToneWin - 1)*10 % 10Hz (40/4 LEDs)
            tbytDat(tt).redLED = []; % omit insufficient number of frames
            tbytDat(tt).redLEDTrainI = [];
            tbytDat(tt).redLEDPulsesOfTrain = [];
        end        
    end

    %% organize CMOS pulses
    % side green CMOS
    if sum(tempSideGreenI)>0
        tbytDat(tt).sideGreenExp = evtInS.sideGreen(tempSideGreenI, 1); % cmos pulses train Id
        tempSideGreenTrainI = unique(evtInS.sideGreen(tempSideGreenI, 2));
        tbytDat(tt).sideGreenExpTrainI = tempSideGreenTrainI;
        tempSideGreenTrainPulses = evtInS.sideGreen(evtInS.sideGreen(:, 2)== tempSideGreenTrainI, 1);
        tempFirstPulseInSideGreenTrain = find(ismember(tempSideGreenTrainPulses, tbytDat(tt).sideGreenExp), 1, 'first');
        tbytDat(tt).sideGreenExpPulsesOfTrain = {tempFirstPulseInSideGreenTrain, tempFirstPulseInSideGreenTrain+numel(tbytDat(tt).sideGreenExp)-1};
        if length(tbytDat(tt).sideGreenExp) < (tbytDat(tt).evtOff-tbytDat(tt).evtOn + p.Results.preToneWin + p.Results.postToneWin - 1)*40 % 40Hz
            tbytDat(tt).sideGreenExp = []; % omit insufficient number of frames
            tbytDat(tt).sideGreenExpTrainI = [];
            tbytDat(tt).sideGreenExpPulsesOfTrain = [];
        end
    end

    % top red CMOS
    if sum(tempTopRedI)>0
        tbytDat(tt).topRedExp = evtInS.topRed(tempTopRedI, 1); % cmos pulses train Id
        tempTopRedTrainI = unique(evtInS.topRed(tempTopRedI, 2));
        tbytDat(tt).topRedExpTrainI = tempTopRedTrainI;
        tempTopRedTrainPulses = evtInS.topRed(evtInS.topRed(:, 2)== tempTopRedTrainI, 1);
        tempFirstPulseInTopRedTrain = find(ismember(tempTopRedTrainPulses, tbytDat(tt).topRedExp), 1, 'first');
        tbytDat(tt).topRedExpPulsesOfTrain = {tempFirstPulseInTopRedTrain, tempFirstPulseInTopRedTrain+numel(tbytDat(tt).topRedExp)-1};
        if length(tbytDat(tt).topRedExp) < (tbytDat(tt).evtOff-tbytDat(tt).evtOn + p.Results.preToneWin + p.Results.postToneWin - 1)*40 % 40Hz
            tbytDat(tt).topRedExp = []; % omit insufficient number of frames
            tbytDat(tt).topRedExpTrainI = [];
            tbytDat(tt).topRedExpPulsesOfTrain = [];
        end
    end
    
    %% organize faceCam pulses
    tbytDat(tt).faceCam = detectFaceCamPulses(tbytDat(tt).evtOn-p.Results.preToneWin, tbytDat(tt).evtOff+p.Results.postToneWin); % faceCam
end

% trial-by-trial peri-stim Lick
tbytLickLogic = arrayfun(@(a, b) a-5 <= evtInS.lick & evtInS.lick <= b+5, [tbytDat(:).evtOn], [tbytDat(:).evtOff], 'un', 0);
%tbytLickLogic = arrayfun(@(a, b) a <= evtInS.lick & evtInS.lick <= b+5, [tbytDat(:).evtOff], [tbytDat(:).evtOff], 'un', 0);
tbytLick = cellfun(@(a) evtInS.lick(a), tbytLickLogic, 'un', 0);
[tbytDat(1:numbTr).Lick] = deal(tbytLick{:});

% water
tbytWaterLogic = arrayfun(@(a, b) a <= evtInS.water & evtInS.water <= b+4, [tbytDat(:).evtOn], [tbytDat(:).evtOff], 'un', 0);
tbytWater = cellfun(@(a) evtInS.water(a), tbytWaterLogic, 'un', 0);
[tbytDat(1:numbTr).water] = deal(tbytWater{:});

% airpuff
tbytAirpuffLogic = arrayfun(@(a, b) a <= evtInS.airpuff & evtInS.airpuff <= b+4, [tbytDat(:).evtOn], [tbytDat(:).evtOff], 'un', 0);
tbytAirpuff = cellfun(@(a) evtInS.airpuff(a), tbytAirpuffLogic, 'un', 0);
[tbytDat(1:numbTr).airpuff] = deal(tbytAirpuff{:});

% trial type info
tt = parse_gng_trialtype(stimopts);
ttName = fieldnames(tt); 
for ff = 1:length(ttName) 
    ttC = num2cell(tt.(ttName{ff})); 
    [tbytDat(1:numbTr).(ttName{ff})] = deal(ttC{:}); 
end

%% Save
if exist(fullfile(filePath, 'Matfiles'), 'dir')==0
    mkdir(fullfile(filePath, 'Matfiles'))
end

header = extract_date_animalID_header(filePath); 

save(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat.mat')), 'tbytDat', 'p');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parseInput_tbyt_auditory_gng_behavior( filePath, vargs )
        % parse input, and extract name-value pairs
        default_preToneWin = 1; % pre-Tone period for event window
        default_postToneWin = 4; % post-Tone window period for event window

        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'preToneWin', default_preToneWin);
        addParameter(p,'postToneWin', default_postToneWin);

        parse(p, filePath, vargs{:})
    end

    function cmosExposureEventsC = organize_cmos_exposure_pulses(cmosExposureEvents, frameShortageAllowance)
        % cmosExposureEvents is assumed to be two columns:
        %   1) time in sec each pulse
        %   2) trial indices

        cmosExposureEventsC = cell(length(unique(cmosExposureEvents(:, 2))), 1);
        for t = 1:length(unique(cmosExposureEvents(:, 2)))
            cmosExposureEventsC{t} = cmosExposureEvents(cmosExposureEvents(:, 2)==t, 1);
        end

        meanFrames = mean(cell2mat(cellfun(@length, cmosExposureEventsC, 'un', 0)));
        enoughFrameLogic = cell2mat(cellfun(@(a) length(a)>=meanFrames-frameShortageAllowance, cmosExposureEventsC, 'un', 0));

        cmosExposureEventsC = cmosExposureEventsC(enoughFrameLogic);

    end

    function header = extract_date_animalID_header(filepath)
        % Extract the portion of the file path that includes the preceding
        % string and the 6-digit date following the underscore.
        % Input:
        %   filepath - the input file path string
        % Output:
        %   header - the extracted string 'preceding_string_6digit_date'

        % Define the regex pattern: look for '\anyString_6digit' pattern
        pattern = '[\\]([a-zA-Z0-9]+_[0-9]{6})[\\]';

        % Use the regexp function to search and return the matched string
        match = regexp(filepath, pattern, 'tokens');

        % If a match is found, output the first match
        if ~isempty(match)
            header = match{1}{1};
        else
            header = ''; % Return an empty string if no match is found
        end
    end
end
