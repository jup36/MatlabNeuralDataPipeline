function parse_tbyt_passive_behavior_without_photodiode(filePath)
% e.g., filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA006/DA006_083023';

%% locate the nidq folder
filePath_nidq = GrabFiles_sort_trials('_g', 0, {filePath});
if isempty(filePath_nidq)
    filePath_nidq = uigetdir(filePath, 'Select the nidq folder');
end

if length(filePath_nidq) > 1
    warning("More than 1 nidq folders found! Will run with the first one!")
end
filePath_nidq = filePath_nidq{1}; % just take the path

%% Load bin file
% binFile = dir(fullfile(filePath_nidq, '*.nidq.bin')); % look for nidq.bin file
% binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
% binFile = binFile(binFileI);
% if length(binFile)>1 || isempty(binFile)
%     error('File could not be found or multiple nidq.bin files exist!');
% end

%% Load evtInS or create one
evtInS = timestamp_behav_events(filePath_nidq, false, 'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'photoDiode'); % behavioral events

%% Load stimopts
filePath_stim = GrabFiles_sort_trials('_stimInfo', 0, {filePath});
if isempty(filePath_stim)
    filePath_stim = uigetdir(filePath, 'Select the stimInfo file!');
end

if length(filePath_stim) > 1
    warning("More than 1 stimInfo files found! Will run with the first one!")
end
filePath_stim = filePath_stim{1}; % just take the path
load(fullfile(filePath_stim), 'stimopts')

%% Identify and sort events
evtInS.water(:, 2) = nan(length(evtInS.water(:, 1)), 1);
evts = evtInS.water;
evts(:, 3) = 3;
assert(size(evts, 1)==sum(stimopts.stim_type==3), 'The number of recorded events and the number of events executed do not match!'); %

% sanity check (match water events before and after sorting)
assert(isequal(evts(evts(:, 3)==3, 1), evtInS.water(:, 1)), 'There is something wrong with event identification and sorting!')

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
[tbytDat(1:numbTr).evtType] = deal(evtTypeC{:}); % 1: common visual, 2: uncommon visual, 3: water reward

% define functions to detect relevant behavioral and task events
detectLicks = @(a, b) evtInS.lick(a <= evtInS.lick & b >= evtInS.lick); % to detect peri-event licks
detectCmosPulses = @(a, b) a <= evtInS.cmosExp(:, 1) & b >= evtInS.cmosExp(:, 1); % to detect cmos exposure pulses
detectFaceCamPulses = @(a, b) evtInS.faceCam(a <= evtInS.faceCam(:, 1) & b >= evtInS.faceCam(:, 1)); % to detect faceCam exposure pulses
detectWater = @(a, b) evtInS.water(a <= evtInS.water(:, 1) & b >= evtInS.water(:, 1)); % to detect water deliveries
%detectAirpuff = @(a, b) evtInS.airpuff(a <= evtInS.airpuff(:, 1) & b >= evtInS.airpuff(:, 1)); % to detect water deliveries

% assign block ID and CMOS exposure pulses
for tt = 1:length(tbytDat)
    % reward trial
    if tbytDat(tt).evtType == 3
        tbytDat(tt).evtOff = tbytDat(tt).evtOn + 4; % 5 seconds were given, take 4 seconds after reward
        % cmos pulses
        tempCmosI = detectCmosPulses(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff); % cmos pulses
        if sum(tempCmosI)>0
            tbytDat(tt).cmosExp = evtInS.cmosExp(tempCmosI, 1); % cmos pulses train Id
            tempCmosTrainI = unique(evtInS.cmosExp(tempCmosI, 2));
            tbytDat(tt).cmosExpTrainI = tempCmosTrainI;
            tempCmosTrainPulses = evtInS.cmosExp(evtInS.cmosExp(:, 2)==tempCmosTrainI, 1);
            tempFirstPulseInTrain = find(ismember(tempCmosTrainPulses, tbytDat(tt).cmosExp), 1, 'first');
            tbytDat(tt).cmosExpPulsesOfTrain = {tempFirstPulseInTrain, tempFirstPulseInTrain+numel(tbytDat(tt).cmosExp)-1};
            if length(tbytDat(tt).cmosExp) < 5*30 % 5s, 30Hz
                tbytDat(tt).cmosExp = []; % omit insufficient number of frames
                tbytDat(tt).cmosExpTrainI = [];
                tbytDat(tt).cmosExpPulsesOfTrain = [];
            end
        end

        % licks
        tbytDat(tt).licks = detectLicks(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff); % licks
        % faceCam
        tbytDat(tt).faceCam = detectFaceCamPulses(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff); % faceCam
        % water
        tbytDat(tt).water = detectWater(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff); % water
        % airpuff
        %tbytDat(tt).airpuff = detectAirpuff(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff); % airpuff

        % visual stimuli trial - note that this session has an issue with photodiode. 
    end
end

%% Save
if exist(fullfile(filePath, 'Matfiles'), 'dir')==0
    mkdir(fullfile(filePath, 'Matfiles'))
end

[~, foldName] = fileparts(filePath);

save(fullfile(filePath, 'Matfiles', strcat(foldName, '_tbytDat.mat')), 'tbytDat');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

end
