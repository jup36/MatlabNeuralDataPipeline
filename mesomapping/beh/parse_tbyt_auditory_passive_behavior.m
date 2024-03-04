function parse_tbyt_auditory_passive_behavior(filePath)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/GC2719/GC2719_020624';

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
evtInS = timestamp_behav_events(filePath_nidq, false, 'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'speaker'); % behavioral events

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
assert(size(evtInS.speaker, 1)==length(stimopts.stim_type), 'The number of recorded events and the number of events executed do not match!')
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
detectTones = @(a, b) a <= evtInS.speaker(:, 1) & b >= evtInS.speaker(:, 1); % to detect tones
detectFaceCamPulses = @(a, b) evtInS.faceCam(a <= evtInS.faceCam(:, 1) & b >= evtInS.faceCam(:, 1)); % to detect faceCam exposure pulses
detectCmosPulses = @(a, b) a <= evtInS.cmosExp(:, 1) & b >= evtInS.cmosExp(:, 1); % to detect cmos exposure pulses
%detectLicks = @(a, b) evtInS.lick(a <= evtInS.lick & b >= evtInS.lick); % to detect peri-event licks
%detectWater = @(a, b) evtInS.water(a <= evtInS.water(:, 1) & b >= evtInS.water(:, 1)); % to detect water deliveries
%detectAirpuff = @(a, b) evtInS.airpuff(a <= evtInS.airpuff(:, 1) & b >= evtInS.airpuff(:, 1)); % to detect water deliveries

% assign block ID and CMOS exposure pulses
for tt = 1:length(tbytDat)
    % auditory stimuli trial
    % cmos pulses
    tempCmosI = detectCmosPulses(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff+1); % cmos pulses
    if sum(tempCmosI)>0
        tbytDat(tt).cmosExp = evtInS.cmosExp(tempCmosI, 1); % cmos pulses train Id
        tempCmosTrainI = unique(evtInS.cmosExp(tempCmosI, 2));
        tbytDat(tt).cmosExpTrainI = tempCmosTrainI;
        tempCmosTrainPulses = evtInS.cmosExp(evtInS.cmosExp(:, 2)==tempCmosTrainI, 1);
        tempFirstPulseInTrain = find(ismember(tempCmosTrainPulses, tbytDat(tt).cmosExp), 1, 'first');
        tbytDat(tt).cmosExpPulsesOfTrain = {tempFirstPulseInTrain, tempFirstPulseInTrain+numel(tbytDat(tt).cmosExp)-1};
        if length(tbytDat(tt).cmosExp) < 4*30 % 4s, 30Hz
            tbytDat(tt).cmosExp = []; % omit insufficient number of frames
            tbytDat(tt).cmosExpTrainI = [];
            tbytDat(tt).cmosExpPulsesOfTrain = [];
        end
    end
    % faceCam
    tbytDat(tt).faceCam = detectFaceCamPulses(tbytDat(tt).evtOn-1, tbytDat(tt).evtOff+1); % faceCam
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
