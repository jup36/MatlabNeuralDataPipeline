function behaviorTimestampsPP(p)
%behaviorTimestamps

%addpath(genpath(''))
%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718';
%p = parse_input_Js(filePath, varargin ); % parse input
%p = parse_input_Js(filePath, {} ); % use this line instead when running line-by-line
% To convert datenum to datetime use, datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum')

cd(p.Results.filePath)

parts = strsplit(p.Results.filePath, filesep);
parent_path = strjoin(parts(1:end-1), filesep);

% load 'stimopts'
stiminfo_file = dir(fullfile([parent_path, filesep, '*_stiminfo.mat'])); 
if isempty(stiminfo_file)
   [stiminfo_file.name, stiminfo_file.folder] = uigetfile(parent_path); 
end
load(fullfile(stiminfo_file.folder, stiminfo_file.name), 'stimopts')

% check the folder whether there's rez file already
if ~isempty(dir(fullfile(p.Results.filePath,'BehVariablesPP.mat'))) % if the BehVariablesPP.mat file already exists in the filePath
    answer = questdlg('BehVariablesPP.mat already exists, Would you like to replace it?','Choice','Replace','Cancel','Cancel');
    switch answer
        case 'Replace'
        case 'Cancel'
            return
    end
end

binFile = dir(fullfile(p.Results.filePath,'*.nidq.bin')); % look for nidq.bin file
binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
binFile = binFile(binFileI);

if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Parse the corresponding metafile
meta  = ReadMeta(binName, p.Results.filePath); % get the meta data (structure)
channels = textscan(meta.acqMnMaXaDw,'%n %n %n %n','Delimiter',',');
nSamp = floor(SampRate(meta));          % sampling rate (default: 25kHz)
totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds

if ~isempty(dir(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'), 'lick', 'trStart', 'reward', 'trEnd', 'camTrig', 'encodeA', 'encodeB', 'laser' ) % if there are gaincorrectedrawtraces already saved, just load them
else
    % Specify the relevant behavioral channel numbers
    cmosExpCh    = channels{1}+p.Results.cmosExpCh;   % ai0
    cmosTrigCh   = channels{1}+p.Results.cmosTrigCh;  % ai1
    speakerCh    = channels{1}+p.Results.speakerCh;   % ai3
    lickCh       = channels{1}+p.Results.lickCh;      % ai4
    lick2Ch      = channels{1}+p.Results.lick2Ch;     % ai6
    bodyCamCh    = channels{1}+p.Results.bodyCamCh;   % ai5
    faceCamCh    = channels{1}+p.Results.faceCamCh;   % ai7
    photoDiodeCh = channels{1}+p.Results.photoDiodeCh; % ai14
    digitCh       = channels{1}+p.Results.digitCh;     % 17 when there are 16 analog inputs recorded

    % preallocate the behavioral data arrays
    cmosExp    = zeros(1,floor(totalTimeSecs*nSamp));
    cmosTrig  = zeros(1,floor(totalTimeSecs*nSamp));
    speaker    = zeros(1,floor(totalTimeSecs*nSamp));
    lick       = zeros(1,floor(totalTimeSecs*nSamp));
    lick2      = zeros(1,floor(totalTimeSecs*nSamp));
    bodyCam    = zeros(1,floor(totalTimeSecs*nSamp));
    faceCam    = zeros(1,floor(totalTimeSecs*nSamp));
    photoDiode = zeros(1,floor(totalTimeSecs*nSamp));
    digit      = zeros(1,floor(totalTimeSecs*nSamp));

    for i = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        start = floor(i * nSamp);
        last  = floor((i+1) * nSamp);

        tempDataArray   = ReadBin(start, nSamp, meta, binName, p.Results.filePath); % read bin data for each second
        tempCmosExp     = tempDataArray(cmosExpCh, :);    
        tempCmosTrig    = tempDataArray(cmosTrigCh, :); 
        tempSpeaker     = tempDataArray(speakerCh, :);    
        tempLick        = tempDataArray(lickCh, :);       
        tempLick2       = tempDataArray(lick2Ch, :);     
        tempBodyCam     = tempDataArray(bodyCamCh, :);    
        tempFaceCam     = tempDataArray(faceCamCh, :);    
        tempPhotoDiode  = tempDataArray(photoDiodeCh, :); 
        tempDigit       = tempDataArray(digitCh, :); 

        cmosExp(1,start+1:last)     = tempCmosExp; % accumulated the decimated data second-by-second
        cmosTrig(1,start+1:last)    = tempCmosTrig;
        speaker(1,start+1:last)     = tempSpeaker;
        lick(1,start+1:last)        = tempLick;
        lick2(1,start+1:last)       = tempLick2;
        bodyCam(1,start+1:last)     = tempBodyCam;
        faceCam(1,start+1:last)     = tempFaceCam;
        photoDiode(1,start+1:last)  = tempPhotoDiode;
        digit(1, start+1:last)      = tempDigit; 
        fprintf('processed %d\n', i+1)
    end
    clearvars i
    
    %digitBin = dec2bin(digit);
    %water = str2num(digitBin(:, 1)); %'str2double' doesn't work here! 
    %airpuff = str2num(digitBin(:, 2)); 
    %figure; hold on; 
    %plot(water); 
    %plot(airpuff); 
    water = digit==4;   % p0.2 -> 2^2
    airpuff = digit==2; % p0.1 -> 2^1

    % Gain correction for analog channnels 
    cmosExp     = GainCorrectNI(cmosExp, 1, meta);      
    cmosTrig    = GainCorrectNI(cmosTrig, 1, meta); 
    speaker     = GainCorrectNI(speaker, 1, meta);   
    lick        = GainCorrectNI(lick, 1, meta);   
    lick2       = GainCorrectNI(lick2, 1, meta);  
    bodyCam     = GainCorrectNI(bodyCam, 1, meta);
    faceCam     = GainCorrectNI(faceCam, 1, meta);
    photoDiode  = GainCorrectNI(photoDiode, 1, meta);
    clearvars temp*
    save('gainCorrectRawTraces', 'meta', 'cmosExp', 'cmosTrig', 'speaker', 'lick', 'lick2', ...
        'bodyCam', 'faceCam', 'photoDiode', 'digit', 'water', 'airpuff')
    fprintf('completed saving raw traces!\n');
end

%% task event detection
if ~isempty(dir(fullfile(p.Results.filePath,'evtIndices.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'evtIndices.mat'),'trStartIdx','trEndIdx','rwdIdx','lickIdx','evtIdx25k','evtIdx1k') % if there are gaincorrectedrawtraces already saved, just load them
else
    
    % detect outcome deliveries
    figure(100); hold on;  
    rwdIdx = detecteventbythreshold(water, nSamp, 50, 'stdFactor', 1, 'plotRez', true, 'chunkPulses', false, 'correctLongPulse', true); % reward (water) deliveries
    fprintf('completed water reward detection!\n');
    punishIdx = detecteventbythreshold(airpuff, nSamp, 50, 'stdFactor', 1, 'plotRez', true, 'chunkPulses', false, 'correctLongPulse', true); % punishment (airpuff) deliveries
    fprintf('completed airpuff punishment detection!\n');
    hold off; 

    [camTrigRiseIdx, camTrigFallIdx, camPulseTrainIdx] = detecteventbythreshold(camTrig, nSamp, 2, 'stdFactor', 1, 'plotRez',false, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', true); % camera trigger

    [trStartIdx,~,~] = detecteventbythreshold(trStart, nSamp, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
    fprintf('completed trial start detection!\n');
    [trEndIdx,~,~]   = detecteventbythreshold(trEnd, nSamp, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial End
    fprintf('completed trial end detection!\n');

    if length(trStartIdx)==length(trEndIdx)
        if ~unique(trEndIdx - trStartIdx>0)
            error('Trial End and Start indices do not make sense!')
        end
    elseif length(trStartIdx)-length(trEndIdx)==1
        if ~unique(trEndIdx - trStartIdx(1,1:length(trEndIdx))>0)
            error('Trial End and Start indices do not make sense!')
        else
            trStartIdx = trStartIdx(1:end-1);
        end
    elseif length(trEndIdx)-length(trStartIdx)==1
        if ~unique(trEndIdx(1,1:length(trStartIdx)) - trStartIdx>0)
            error('Trial End and Start indices do not make sense!')
        else
            trEndIdx = trEndIdx(1:end-1);
        end
    else
        error('Trial End and Start indices do not make sense!')
    end

    rwdIdx  = detecteventbythreshold(reward, nSamp, 50, 'stdFactor',1, 'plotRez',false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse',true);  % reward
    fprintf('completed reward detection!\n');

    % detect licks
    %lickIdx = detecteventbythreshold(decimate(lick,nSamp/1000), 1000, 20, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false); % lick, detect licks after downsampling (there seems to be some denoising effect with decimation, see /Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718/lickSignalExample_1kHzVS25kHz.fig as an example)
    [lickIdx_decimated, ~] = detecteventbythreshold_with_artifact_correction(lick, nSamp, true);
    fprintf('completed lick detection!\n');

    [camTrigRiseIdx, camTrigFallIdx, camPulseTrainIdx] = detecteventbythreshold(camTrig, nSamp, 2, 'stdFactor', 1, 'plotRez',false, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', true); % camera trigger
    %camTrigRiseIdx1ms = round(camTrigRiseIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
    %camTrigFallIdx1ms = round(camTrigFallIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
    fprintf('completed camera trigger detection!\n');

    % store evt indices
    evtIdx25k.trStartIdx = trStartIdx;
    evtIdx25k.trEndIdx = trEndIdx;
    evtIdx25k.rwdIdx  = rwdIdx+nSamp*(p.Results.rewardDelay/1000); % correct for the delay
    evtIdx25k.lickIdx = 'NA';
    evtIdx25k.camTrigRiseIdx = camTrigRiseIdx;
    evtIdx25k.camTrigFallIdx = camTrigFallIdx;
    evtIdx25k.camPulseTrainIdx = camPulseTrainIdx;

    evtIdx1k.trStartIdx = round(trStartIdx./25);
    evtIdx1k.trEndIdx = round(trEndIdx./25);
    evtIdx1k.rwdIdx  = round(rwdIdx./25)+p.Results.rewardDelay;
    evtIdx1k.lickIdx = round(lickIdx_decimated); % no need to divide it again if already decimated
    evtIdx1k.camTrigRiseIdx = round(camTrigRiseIdx./25);
    evtIdx1k.camTrigFallIdx = round(camTrigFallIdx./25);
    evtIdx1k.camPulseTrainIdx = camPulseTrainIdx; % pulse Train Id

    save('evtIndices','trStartIdx','trEndIdx','rwdIdx','lickIdx_decimated','evtIdx25k','evtIdx1k')
end

%% spot the trial-by-trial and all-trials behavioral csv files
if isempty(dir(fullfile(p.Results.filePath,'20*')))
    error('Cannot find the trial-by-trial behavior data csv files!')
end

behFilePath = dir(fullfile(p.Results.filePath,'20*')); % dir where the trial-by-trial behavioral csv files are saved
%tbytCsvList = dir(fullfile(behFilePath.folder,behFilePath.name,'trial_*'));    % trial-by-trial files
allTrialCsv = dir(fullfile(behFilePath.folder,behFilePath.name,'trials.csv')); % all trial file
if length(allTrialCsv)==1
    trialsFileName = fullfile(allTrialCsv.folder,allTrialCsv.name);
    trialsCsv = readtable(trialsFileName);
else
    error('More than one trials.csv file detected!')
end


if size(trialsCsv,1) < length(trStartIdx)
    if length(trStartIdx) - size(trialsCsv,1) == 1 && length(trEndIdx) - size(trialsCsv,1) == 1
        trStartIdx = trStartIdx(1:size(trialsCsv,1));
        trEndIdx = trEndIdx(1:size(trialsCsv,1));
    else
        error('The # of trials does not match between nidq and csv!')
    end
end

%% parse stepper encoder data; pin A, pin B, and extract joystick kinematics
jsDistCoeff = 2*pi*90/4000; % joystick movement distance conversion coefficient (4000: the number of total edges(rises/falls) per resolution of the encoder)
% get the midpoint of the max and min of step and direc as the mean of 10 folds to deal with outliers
tenFold = linspace(1,length(encodeA),10);
pinAmaxFolds = []; pinAminFolds = []; % pinA max and min of 10 folds
pinBmaxFolds = []; pinBminFolds = []; % pinB max of 10 folds
for i=1:10 % 10 folds
    if i<10
        pinAmaxFolds = [pinAmaxFolds, nanmax(encodeA(1,floor(tenFold(i)):floor(tenFold(i+1))-1))];
        pinAminFolds = [pinAminFolds, nanmin(encodeA(1,floor(tenFold(i)):floor(tenFold(i+1))-1))];
        pinBmaxFolds = [pinBmaxFolds, nanmax(encodeB(1,floor(tenFold(i)):floor(tenFold(i+1))-1))];
        pinBminFolds = [pinBminFolds, nanmin(encodeB(1,floor(tenFold(i)):floor(tenFold(i+1))-1))];
    elseif i==10
        pinAmidPoint = (nanmean(pinAmaxFolds)+nanmean(pinAminFolds))/2;
        pinBmidPoint = (nanmean(pinBmaxFolds)+nanmean(pinBminFolds))/2;
    end
end
clearvars i

jsTime25k = struct;
jsTime1k = struct;

% get the trial information reconstructed based off of trialsCsv, trStartIdx, and rwdIdx
% [trialInfo] = reconstructAssay( trialsCsv, trStartIdx, rwdIdx );
trialInfo.reachPos1 = trialsCsv.block_reach_position_1;
trialInfo.pullTq = trialsCsv.block_pull_torque;
trialInfo.pullThreshold = trialsCsv.pull_threshold;


for t = 1:length(trStartIdx) % increment trials
    if ~isempty(find(trEndIdx>trStartIdx(t),1)) % if there's a trEnd
        jsTime25k(t).trStart = trStartIdx(t); % trStart detected by the go-cue onset
        % redefine the trial start as the joystick in position timepoint (trJsReadyTime) by examining the baseline encoder data,
        % as the trStart defined with the cue tone onset doesn't align perfectly with the joystick being in the start position.
        jsTime25k(t).trEnd = trEndIdx(find(trEndIdx > trStartIdx(t),1,'first'));  % get the trEnd time point

        trBaseRange = trStartIdx(t)-2*nSamp:trStartIdx(t)+(nSamp/10/2)-1; % 2 sec baseline, just take 2 sec before the cue onset, with padding on the righthand side corresponding to 50 ms
        trBasePinA = encodeA(trBaseRange); % pin A this trial baseline
        trBasePinB = encodeB(trBaseRange); % pin B this trial baseline
        biTrBasePinA = double(trBasePinA>pinAmidPoint); % binarize the pinA signal
        biTrBasePinB = double(trBasePinB>pinBmidPoint); % binarize the pinB signal

        [trBaseJsState, ~] = readStepperEncoder(biTrBasePinA, biTrBasePinB); % read out the Js position from the binarized quadrature encoded signals
        [trBaseLatestStillTime, jsTime1k(t).baseJsTrajmm, jsTime1k(t).baseSmJsVel, jsTime1k(t).basePeriodicAbsVelSum] = findJsReadyPt2(cumsum(trBaseJsState), nSamp); % find the Js ready time by detecting a still point in the later portion of the baseline period

        jsTime25k(t).pull_threshold = trialInfo.pullThreshold(t); % pull threshold
        jsTime25k(t).pull_torque = trialInfo.pullTq(t); % pull torque
        jsTime25k(t).reachP1 = trialInfo.reachPos1(t); % reach position 1

        if isnan(trBaseLatestStillTime)
            jsTime25k(t).trJsReady = nan;
            jsTime25k(t).trJsTraj = nan;
            jsTime25k(t).dctrJsTraj = nan;
            jsTime25k(t).rewardT = nan;
            % determine if this trial led to reward
            if ~isempty(find(abs(rwdIdx-jsTime25k(t).trEnd)<=nSamp,1)) % in case there's reward delivery within 1-sec window relative to the trial end
                jsTime25k(t).rewarded = true;
            else
                jsTime25k(t).rewarded = false;
            end

            jsTime25k(t).trialType = nan; % trialType
            jsTime1k(t).movKins = nan; % movement kinematics assign nan

        elseif trStartIdx(t)-2*nSamp+trBaseLatestStillTime <= jsTime25k(t).trStart+(nSamp/10/2) % if there's a joystick still point within 50 ms from the trStart signal
            jsTime25k(t).trJsReady = trStartIdx(t)-2*nSamp+trBaseLatestStillTime; % get the trJsReady time

            trRange = jsTime25k(t).trJsReady:jsTime25k(t).trEnd; % trial range aligned to the Js ready time

            trPinA = encodeA(trRange); % pin A this trial aligned to the Js ready
            trPinB = encodeB(trRange); % pin A this trial aligned to the Js ready

            biTrPinA = double(trPinA>pinAmidPoint); % binarize the encoder signal
            biTrPinB = double(trPinB>pinBmidPoint); % binarize the encoder signal
            [trJsState, ~] = readStepperEncoder(biTrPinA, biTrPinB); % read out the Js position from the binarized quadrature encoded signals
            jsTime25k(t).trJsTraj = cumsum(trJsState);  % get the cumulative joystick trajectory for this trial
            jsTime25k(t).dctrJsTraj = decimate(jsTime25k(t).trJsTraj ,round(nSamp/1000)); % decimate the Traj (downsampling the 25kHz data at 1kHz, so the time resoluation to be 1ms)

            % to use sg filter get sgfiltFramelen
            if length(jsTime25k(t).dctrJsTraj) >= p.Results.sgfiltFramelen
                sgfiltFramelen = p.Results.sgfiltFramelen;
            elseif length(jsTime25k(t).dctrJsTraj) < p.Results.sgfiltFramelen
                if mod(length(jsTime25k(t).dctrJsTraj),2)==0
                    sgfiltFramelen = length(jsTime25k(t).dctrJsTraj)-1; % the frame length for sg filter needs to be an odd number
                else
                    sgfiltFramelen = length(jsTime25k(t).dctrJsTraj);
                end
            end

            tempSmdctrJsTraj = sgolayfilt(jsTime25k(t).dctrJsTraj,3,sgfiltFramelen);
            %tempSmdctrJsTrajmm = tempSmdctrJsTraj*jsDistCoeff; % convert the Js traj into mm

            % determine if the trial got rewarded
            if ~isempty(find(abs(rwdIdx-jsTime25k(t).trEnd)<=nSamp,1)) % in case there's reward delivery within 1-sec window relative to the trial end
                jsTime25k(t).rewarded = true;
            else
                jsTime25k(t).rewarded = false;
            end

            % classify the trial ('sp': successfull pull, 'ps': push, 'pm': premature pull, 'to': timeout, 'nn': not identified)
            if jsTime25k(t).rewarded % if rewarded
                if ~isempty(find(tempSmdctrJsTraj<jsTime25k(t).pull_threshold,1)) % check the negative threshold crossing
                    jsTime25k(t).trialType = 'sp'; % successful pull
                    jsTime25k(t).rewardT = jsTime25k(t).trEnd+nSamp*(p.Results.rewardDelay/1000); % reward if delivered must be 1000ms (or some other delay) after the trial end
                else % this should be an premature pull which has been accidentally rewarded
                    jsTime25k(t).trialType = nan; % just throw out these trials for now
                    jsTime25k(t).rewardT = nan;
                end
            else % if not rewarded
                jsTime25k(t).rewardT = nan;
                if length(jsTime25k(t).dctrJsTraj)>p.Results.trialTimeout-100 % timeout
                    jsTime25k(t).trialType = 'to';
                elseif isempty(find(jsTime25k(t).dctrJsTraj<jsTime25k(t).pull_threshold,1)) && ~isempty(find(jsTime25k(t).dctrJsTraj>50,1)) % push (no pull beyond the pull threshold && push beyond a certain threshold)
                    jsTime25k(t).trialType = 'ps';
                elseif ~isempty(find(jsTime25k(t).dctrJsTraj<jsTime25k(t).pull_threshold,1)) % pull (unrewarded)
                    impullThresCross = find(jsTime25k(t).dctrJsTraj<jsTime25k(t).pull_threshold,1); % premature pull threshold crossing point
                    if ~isempty(find(jsTime25k(t).dctrJsTraj(impullThresCross:end)>0,1))
                        jsTime25k(t).trialType = 'pmpp'; % premature pull and push
                    else
                        jsTime25k(t).trialType = 'pm'; % premature pull
                    end
                else
                    jsTime25k(t).trialType = 'nn';
                end
            end

            % trial-by-trial reach kinematics analysis
            if unique(jsTime25k(t).pull_torque==p.Results.meanMass(1,:))==false % try interpolation
                inter_meanMass(1, :) = p.Results.meanMass(1, 1):1:p.Results.meanMass(1, end);
                inter_meanMass(2, :) = interp1(p.Results.meanMass(1, :), p.Results.meanMass(2, :), p.Results.meanMass(1, 1):1:p.Results.meanMass(1, end));
                tempMass = inter_meanMass(2,jsTime25k(t).pull_torque==inter_meanMass(1, :));
            else % no need to interpolate
                tempMass = p.Results.meanMass(2,jsTime25k(t).pull_torque==p.Results.meanMass(1,:));
            end
            [ jsTime1k(t).movKins ] = jsReachKinematics( jsTime25k(t).dctrJsTraj, jsTime25k(t).pull_threshold, jsTime25k(t).trialType, tempMass, sgfiltFramelen );
            %figure; plot(interp1(1:length(jsTime25k(t).csvTraj),jsTime25k(t).csvTraj,1:length(jsTime25k(t).dcsmtrJsTraj)));
        end
    else % if there's no trialEnd left
    end
    fprintf('completed trial #%d\n', t);
end
clearvars temp* t

% transfer relevant information to jsTime1k struct from jsTime25k
n2cTrStart = num2cell(round([jsTime25k(:).trStart]./25)); [jsTime1k.trStart] = n2cTrStart{:};
n2cTrEnd   = num2cell(round([jsTime25k(:).trEnd]./25));   [jsTime1k.trEnd]   = n2cTrEnd{:};
n2cTrJsReady = num2cell(round([jsTime25k(:).trJsReady]./25)); [jsTime1k.trJsReady] = n2cTrJsReady{:};
n2cRewarded = num2cell([jsTime25k(:).rewarded]); [jsTime1k.rewarded] = n2cRewarded{:};
n2cPull_threshold = num2cell([jsTime25k(:).pull_threshold]); [jsTime1k.pull_threshold] = n2cPull_threshold{:};
n2cPull_torque = num2cell([jsTime25k(:).pull_torque]); [jsTime1k.pull_torque] = n2cPull_torque{:};
n2cReachP1 = num2cell([jsTime25k(:).reachP1]); [jsTime1k.reachP1] = n2cReachP1{:};
n2cRewardT = num2cell(round([jsTime25k(:).rewardT]./25)); [jsTime1k.rewardT] = n2cRewardT{:};
if p.Results.laserUsed
    n2cStimLaserOn = num2cell(round([jsTime25k(:).stimLaserOn]./25)); [jsTime1k.stimLaserOn] = n2cStimLaserOn{:};
    n2cStimLaserOff = num2cell(round([jsTime25k(:).stimLaserOff]./25)); [jsTime1k.stimLaserOff] = n2cStimLaserOff{:};
end
[jsTime1k.trialType] = jsTime25k(:).trialType;

% generate a plot to inspect a certain trial
%movKinsPlot(jsTime1k(t).movKins);

%% Save relevant BehVariables
cd(p.Results.filePath)
save('BehVariablesPP', 'jsTime1k', 'jsTime25k', 'evtIdx25k', 'evtIdx1k', 'p', 'trialsCsv', 'trialInfo'); % 'tbytCsvList') % append the position/velocity data variables

end
