function evtInS = timestamp_behav_events(path_raw, redetect_logic, varargin)
%This function takes the number of behavioral variables as varargin and
% gets the timestamps (in sec) of the variables and saves them as a
% structure named evtInS in the files origianl nidq directory.
% TO DO: Add more variables in the SWITCH CASE block.

% varargin = {'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'photoDiode', 'speaker'}; 

% get stimInfo
filePath = fileparts(path_raw); 
filePath_stimInfo = GrabFiles_sort_trials('_stimInfo', 0, {filePath});
if isempty(filePath_stimInfo)
    filePath_stimInfo = GrabFiles_sort_trials('_stimInfo', 0, {fullfile(filePath, 'Matfiles')});
    if isempty(filePath_stimInfo)
         [stimInfo_file, stimInfo_folder] = uigetfile('*_stimInfo*', 'Select the stimInfo file', filePath);
         filePath_stimInfo = {fullfile(stimInfo_folder, stimInfo_file)};    
    end
end
load(filePath_stimInfo{1}, 'stimopts'); 

% Load bin file
binFile = dir(fullfile(path_raw, '*.nidq.bin')); % look for nidq.bin file
binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
binFile = binFile(binFileI);

if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Load evtInS
if redetect_logic
    evtInS = struct;
else
    if ~isempty(findFileWithString(path_raw, 'evtInS'))
        load(fullfile(path_raw, 'evtInS'))
        evtInS_fields = fieldnames(evtInS);
        evtInS_isfield = cell2mat(cellfun(@(a) any(strcmp(a, evtInS_fields)), varargin, 'un', 0));
        if sum(evtInS_isfield)==length(varargin)
            return; % if all the events were already detected and available just return
        else % in case there's still to be detected and appended
            varargin = varargin(~evtInS_isfield);
        end
    else
        evtInS = struct;
    end
end

% Parse the corresponding metafile
meta  = ReadMeta(binName, path_raw); % get the meta data (structure)
sample_rate=str2double(meta.niSampRate);
nSamp = floor(SampRate(meta));          % sampling rate (default: 25kHz)

output = struct;

%% Main (load raw traces and perform event detection)
if exist(fullfile(path_raw,'gainCorrectRawTraces.mat'), 'file')~=2 || redetect_logic
    p = parse_input_PP(path_raw, {});
    behaviorTimestampsPP(p);
end

if exist(fullfile(path_raw,'gainCorrectRawTraces.mat'), 'file')==2
    for jj = 1:length(varargin)
        fieldName = sprintf('variable_%d', jj);
        output.(fieldName) = load(fullfile(path_raw,'gainCorrectRawTraces.mat'), varargin{jj});
        assert(isfield(output.(fieldName), varargin{jj})) % ensure the raw trace is loaded properly

        switch varargin{jj}
            case 'lick'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'lick'))
                    lick = output.(fieldName).lick;
                    timeStamp = (1:length(lick))';
                    timeStampSec = timeStamp./sample_rate;
                    lickIdx = detecteventbythreshold_noAbs(path_raw, lick, nSamp, 50, 'lick', ...
                        'stdFactor', 3, 'plotRez', true, 'chunkPulses', false, 'correctLongPulse', true);
                    evtInS.lick = timeStampSec(lickIdx);
                    fprintf('completed lick detection!\n');
                end
            case 'faceCam'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'faceCam'))
                    faceCam = output.(fieldName).faceCam;
                    timeStamp = (1:length(faceCam))';
                    timeStampSec = timeStamp./sample_rate;
                    [faceCamRiseIdx, ~, faceCamTrainIdx] = detecteventbythreshold_interactive(path_raw, faceCam, nSamp, 4, 'faceCam', ...
                        'stdFactor', 2, 'plotRez', true, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', false, ...
                        'filterArtifact', true, 'lowpassCutoff', 5, 'nyquist', 100, 'removeOffTrainPulses', true, ...
                        'detectLater', 5, 'detectEarlier', 5); % camera trigger, detect only after 5s also cutting off the final 5 s. 
                    evtInS.faceCam = [timeStampSec(faceCamRiseIdx), faceCamTrainIdx'];
                    fprintf('completed faceCam detection!\n');
                end
            case 'cmosExp'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'cmosExp'))
                    cmosExp = output.(fieldName).cmosExp;
                    timeStamp = (1:length(cmosExp))';
                    timeStampSec = timeStamp./sample_rate;
                    [cmosExpRiseIdx, ~, cmosExpTrainIdx] = detecteventbythreshold_interactive(path_raw, cmosExp, nSamp, 20, 'cmosExp', ...
                        'stdFactor', 1, 'plotRez', true,... 
                        'chunkPulses', true, 'chunkInterval', 500, 'correctLongPulse', false, ...
                        'filterArtifact', true, 'lowpassCutoff', 1, 'nyquist', 100, 'removeOffTrainPulses', true, ...
                        'detectLater', 5, 'detectEarlier', 5); % CMOS trigger
                    evtInS.cmosExp = [timeStampSec(cmosExpRiseIdx), cmosExpTrainIdx'];
                    fprintf('completed cmosExp detection!\n');
                end
            case 'water'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'water'))
                    water = output.(fieldName).water;
                    timeStamp = (1:length(water))';
                    timeStampSec = timeStamp./sample_rate;
                    waterIdx = detecteventbythreshold_noAbs(path_raw, water, nSamp, 50, 'water', ...
                        'stdFactor', 1.5, 'plotRez', true, 'chunkPulses', false, 'correctLongPulse', true); % CMOS trigger
                    evtInS.water = timeStampSec(waterIdx);
                    fprintf('completed water detection!\n');
                end
            case 'photoDiode'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'photoDiode'))
                    photoDiode = output.(fieldName).photoDiode;
                    timeStamp = (1:length(photoDiode))';
                    timeStampSec = timeStamp./sample_rate;
                    [photoDiodeRiseIdx, photoDiodeFallIdx ] = detecteventbythreshold_interactive_timeout(path_raw, photoDiode, nSamp, 50, 'photoDiode', ...
                        'stdFactor', 1, 'plotRez', true,... 
                        'chunkPulses', false, 'correctLongPulse', false, 'findEarlyOnset', true, ...
                        'filterArtifact', false, 'earlyCutoff', 1, 'lateCutoff', length(photoDiode)); % CMOS trigger
                    evtInS.photoDiode = [timeStampSec(photoDiodeRiseIdx), timeStampSec(photoDiodeFallIdx)];
                    fprintf('completed photoDiode detection!\n');
                end
            case 'airpuff'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'airpuff'))
                    airpuff = output.(fieldName).airpuff;
                    if isempty(airpuff)
                        evtInS.airpuff = []; 
                    else
                        timeStamp = (1:length(airpuff))';
                        timeStampSec = timeStamp./sample_rate;
                        airpuffIdx = detecteventbythreshold_noAbs(path_raw, airpuff, nSamp, 50, 'airpuff', ...
                            'stdFactor', 1.5, 'plotRez', true, 'chunkPulses', false, 'correctLongPulse', true); % CMOS trigger
                        evtInS.airpuff = timeStampSec(airpuffIdx);
                        fprintf('completed airpuff detection!\n');
                    end
                end
            case 'speaker'
                if ~(exist('evtInS', 'var')==1 && isfield(evtInS, 'speaker'))
                    speaker = output.(fieldName).speaker;
                    timeStamp = (1:length(speaker))';
                    timeStampSec = timeStamp./sample_rate;
                    [speakerRiseIdx, speakerFallIdx ] = detecteventbythreshold_interactive_timeout(path_raw, speaker, nSamp, 2000, 'speaker', ...
                        'stdFactor', 1, 'plotRez', true, ... 
                        'chunkPulses', false, 'correctLongPulse', false, 'findEarlyOnset', false, ...
                        'filterArtifact', false, 'earlyCutoff', 1, 'lateCutoff', length(speaker), ...
                        'useAudioConvolution', true, 'freqForConvolution', stimopts.tone_table(:, 1)); % speaker trigger detect using convoluted signals
                    evtInS.speaker = [timeStampSec(speakerRiseIdx), timeStampSec(speakerFallIdx)];
                    fprintf('completed tone detection!\n');
                end
        end
    end
end

%% save evtInS
if exist(fullfile(path_raw, "evtInS"), 'file') == 2
    save(fullfile(path_raw, "evtInS"), 'evtInS', '-append')
else
    save(fullfile(path_raw, "evtInS"), "evtInS")
end

fprintf('completed event detection and saved events in the nidq folder!\n');

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = parse_input_PP( filePath, vargs )
% parse input, and extract name-value pairs
default_reReadBin = false; % by default do not re-read the raw bin file, if done already
default_numbNeuralProbe = 0;  % specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
default_numbChEachProbe = 64; % specify how many channels are on the probe
default_cmosExpCh   = 1;    % ai0: CMOS exposure pulses
default_cmosTrigCh  = 2;    % ai1: CMOS trigger
default_speakerCh   = 4;    % ai3: reward delivery
default_lickCh      = 5;    % ai4: lick detector1
default_lick2Ch     = 7;    % ai6: lick detector2
default_bodyCamCh   = 6;    % ai5: bodyCam
default_faceCamCh   = 8;    % ai7: faceCam
default_photoDiodeCh = 15;  % ai14: photo diode
default_digitCh     = 17;   % the digital channel is always at the tailend

p = inputParser; % create parser object
addRequired(p,'filePath');
addParameter(p,'reReadBin',default_reReadBin);
addParameter(p,'numbNeuralProbe',default_numbNeuralProbe);
addParameter(p,'numbChEachProbe',default_numbChEachProbe);
addParameter(p,'cmosExpCh',default_cmosExpCh);
addParameter(p,'cmosTrigCh',default_cmosTrigCh);
addParameter(p,'speakerCh',default_speakerCh);
addParameter(p,'lickCh',default_lickCh);
addParameter(p,'lick2Ch',default_lick2Ch);
addParameter(p,'bodyCamCh',default_bodyCamCh);
addParameter(p,'faceCamCh',default_faceCamCh);
addParameter(p,'photoDiodeCh',default_photoDiodeCh);
addParameter(p,'digitCh',default_digitCh);

parse(p,filePath,vargs{:})
end