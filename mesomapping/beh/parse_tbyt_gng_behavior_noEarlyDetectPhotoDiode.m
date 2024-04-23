function parse_tbyt_gng_behavior_noEarlyDetectPhotoDiode(filePath)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_101523';

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
%evtInS = timestamp_behav_events(filePath_nidq, true, 'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'photoDiode'); % behavioral events
evtInS = timestamp_behav_events_noEarlyDetectPhotoDiode(filePath_nidq, true,'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'photoDiode'); % behavioral events

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
assert(size(evtInS.photoDiode, 1)==length(stimopts.stim_type), 'The number of recorded events and the number of events executed do not match!')

%% tbytDat
numbTr = size(evtInS.photoDiode, 1);
tbytDat = struct;
% stim On
stimOnC = num2cell(evtInS.photoDiode(:, 1));
[tbytDat(1:numbTr).stimOn] = deal(stimOnC{:});
% stim Off
stimOffC = num2cell(evtInS.photoDiode(:, 2));
[tbytDat(1:numbTr).stimOff] = deal(stimOffC{:});

% trial-by-trial peri-stim Lick
tbytLickLogic = arrayfun(@(a, b) a-5 <= evtInS.lick & evtInS.lick <= b+5, [tbytDat(:).stimOn], [tbytDat(:).stimOff], 'un', 0);
tbytLick = cellfun(@(a) evtInS.lick(a), tbytLickLogic, 'un', 0);
[tbytDat(1:numbTr).Lick] = deal(tbytLick{:});

% trial-by-trial during-stim Lick
tbytStimLickLogic = arrayfun(@(a, b) a <= evtInS.lick & evtInS.lick <= b, [tbytDat(:).stimOn], [tbytDat(:).stimOff], 'un', 0);
tbytStimLick = cellfun(@(a) evtInS.lick(a), tbytStimLickLogic, 'un', 0);
[tbytDat(1:numbTr).stimLick] = deal(tbytStimLick{:});

% water
tbytWaterLogic = arrayfun(@(a, b) a <= evtInS.water & evtInS.water <= b+2, [tbytDat(:).stimOn], [tbytDat(:).stimOff], 'un', 0);
tbytWater = cellfun(@(a) evtInS.water(a), tbytWaterLogic, 'un', 0);
[tbytDat(1:numbTr).water] = deal(tbytWater{:});

% airpuff
tbytAirpuffLogic = arrayfun(@(a, b) a <= evtInS.airpuff & evtInS.airpuff <= b+2, [tbytDat(:).stimOn], [tbytDat(:).stimOff], 'un', 0);
tbytAirpuff = cellfun(@(a) evtInS.airpuff(a), tbytAirpuffLogic, 'un', 0);
[tbytDat(1:numbTr).airpuff] = deal(tbytAirpuff{:});

% CMOS exposure
% cmosEvtC = organize_cmos_exposure_pulses(evtInS.cmosExp, 10);
% if numbTr == numel(cmosEvtC)
%     [tbytDat(1:numbTr).cmos] = deal(cmosEvtC{:});
% else
%     warning('The # of cmos exposure events does not match with the # of trials!')
% end

% faceCam exposure
faceCamEvtC = organize_cmos_exposure_pulses(evtInS.faceCam, 10);
if numbTr == numel(faceCamEvtC)
    [tbytDat(1:numbTr).faceCam] = deal(faceCamEvtC{:});
else
    warning('The # of faceCam exposure events does not match with the # of trials!')
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