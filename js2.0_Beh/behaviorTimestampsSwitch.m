function behaviorTimestampsSwitch(p)
%behaviorTimestamps

%addpath(genpath(''))
%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718';
%p = parse_input_Js(filePath, varargin ); % parse input
%p = parse_input_Js(filePath, {} ); % use this line instead when running line-by-line
% To convert datenum to datetime use, datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum')

cd(p.Results.filePath)

% check the folder whether there's rez file already
if ~isempty(dir(fullfile(p.Results.filePath,'BehVariablesSwitch.mat'))) % if the BehVariablesSwitch.mat file already exists in the filePath
    answer = questdlg('BehVariablesSwitch.mat already exists, Would you like to replace it?','Choice','Replace','Cancel','Cancel');
    switch answer
        case 'Replace'
        case 'Cancel'
            return
    end
end

binFile = dir(fullfile(p.Results.filePath,'*.nidq.bin')); % look for nidq.bin file

if length(binFile)>1 || isempty(binFile)
    disp('Select the directory where the nidq.bin file exists!!')
    binFileDir = uigetdir(p.Results.filePath); 
    if ~isempty(binFileDir) && size(binFileDir,1)==1
        binFile = dir(fullfile(binFileDir,'*.nidq.bin')); % look for nidq.bin file
    else
        error('File could not be found or multiple nidq.bin files exist!');
    end
end

binName = binFile.name;

% Parse the corresponding metafile
meta  = ReadMeta(binName, binFile.folder); % get the meta data (structure)
channels = textscan(meta.acqMnMaXaDw,'%n %n %n %n','Delimiter',',');
nSamp = round(SampRate(meta),0); % sampling rate (default: 8kHz)
%totalTimeSecsMeta = str2double(meta.fileTimeSecs); % total duration of file in seconds

if ~isempty(dir(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'), 'trStart', 'reward', 'trEnd', 'camTrig', 'encodeA', 'encodeB', 'position1', 'position2', 'position3', 'laser' ) % if there are gaincorrectedrawtraces already saved, just load them
else
    % Specify the relevant behavioral channel numbers
    camTrigCh  = channels{1}+p.Results.camTrigCh; % ch# for camera trigger
    trStartCh  = channels{1}+p.Results.trStartCh; % ch# for trial start
    encodeACh  = channels{1}+p.Results.encodeACh; % ch# for stepper direction
    encodeBCh  = channels{1}+p.Results.encodeBCh; % ch# for stepper steps    
    rewardCh   = channels{1}+p.Results.rewardCh; % ch# for reward delivery
    trEndCh    = channels{1}+p.Results.trEndCh;  % ch# for trial end (either by successful pull or error/timeout)
    lickCh     = channels{1}+p.Results.lickCh; % ch# for lick detect
    position1Ch = channels{1}+p.Results.position1Ch; % ch# for position1Ch
    position2Ch = channels{1}+p.Results.position2Ch; % ch# for position2Ch
    laserCh   = channels{1}+p.Results.laserCh;   % ch# for laser triggers
    pacerCh = channels{1}+p.Results.pacerCh;   % ch# for laser triggers
    % file size check
    data_info = dir(fullfile(binFile.folder, binFile.name));
    data_size = data_info.bytes;
    sample_number = data_size / (2 * str2double(meta.nSavedChans));
    totalTimeSecs = floor(sample_number/str2double(meta.niSampRate)); 
    if ~isequal(str2double(meta.fileSizeBytes),data_size)
        warning('Working on a modified bin file!')
        meta.fileSizeBytes = num2str(data_size); 
    end
    
    % preallocate the behavioral data arrays
    camTrig   = zeros(1,floor(totalTimeSecs*8000));
    trStart   = zeros(1,floor(totalTimeSecs*8000));
    encodeA   = zeros(1,floor(totalTimeSecs*8000));
    encodeB   = zeros(1,floor(totalTimeSecs*8000));
    reward    = zeros(1,floor(totalTimeSecs*8000));
    trEnd     = zeros(1,floor(totalTimeSecs*8000));
    lick      = zeros(1,floor(totalTimeSecs*8000));
    position1 = zeros(1,floor(totalTimeSecs*8000));
    position2 = zeros(1,floor(totalTimeSecs*8000));
    laser     = zeros(1,floor(totalTimeSecs*8000));
    pacer     = zeros(1,floor(totalTimeSecs*8000));
    
    for k = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        tempDataArray = ReadBin(k*nSamp, nSamp, meta, binName, binFile.folder); % read bin data for each second
        
        tempCamTrig  = tempDataArray(camTrigCh,:); % do not decimate for higher temporal resolution
        tempTrStart  = tempDataArray(trStartCh,:); % decimate the data
        tempEncodeA  = tempDataArray(encodeACh,:); % do not decimate for higher temporal resolution
        tempEncodeB  = tempDataArray(encodeBCh,:); % do not decimate for higher temporal resolution        
        tempReward   = tempDataArray(rewardCh,:);
        tempTrEnd    = tempDataArray(trEndCh,:);
        tempLick     = tempDataArray(lickCh,:);
        tempPosition1 = tempDataArray(position1Ch,:);
        tempPosition2 = tempDataArray(position2Ch,:);   
        tempLaser    = tempDataArray(laserCh,:);
        tempPacer = tempDataArray(pacerCh,:);
        
        camTrig(1,k*8000+1:(k+1)*8000) = tempCamTrig;
        trStart(1,k*8000+1:(k+1)*8000) = tempTrStart; % accumulated the decimated data second-by-second
        encodeA(1,k*8000+1:(k+1)*8000) = tempEncodeA;
        encodeB(1,k*8000+1:(k+1)*8000) = tempEncodeB;
        reward(1,k*8000+1:(k+1)*8000) = tempReward;
        trEnd(1,k*8000+1:(k+1)*8000)  = tempTrEnd;
        lick(1,k*8000+1:(k+1)*8000) = tempLick;
        position1(1,k*8000+1:(k+1)*8000) = tempPosition1; 
        position2(1,k*8000+1:(k+1)*8000) = tempPosition2; 
        laser(1,k*8000+1:(k+1)*8000) = tempLaser;    
        pacer(1,k*8000+1:(k+1)*8000) = tempPacer;   
        fprintf('processed %d\n', k+1)
    end
    clearvars k
    
    % Gain correction for channnels of interest
    camTrig = GainCorrectNI(camTrig, 1, meta); % gain-corrected voltage trace for camTrig
    trStart = GainCorrectNI(trStart, 1, meta); % gain-corrected voltage trace for trStart
    encodeA = GainCorrectNI(encodeA, 1, meta); % gain-corrected voltage trace for encodeA
    encodeB = GainCorrectNI(encodeB, 1, meta); % gain-corrected voltage trace for encodeB
    reward = GainCorrectNI(reward, 1, meta);   % gain-corrected voltage trace for reward
    trEnd  = GainCorrectNI(trEnd, 1, meta);    % gain-corrected voltage trace for trEnd
    lick = GainCorrectNI(lick, 1, meta); % gain-corrected voltage trace for lick
    position1 = GainCorrectNI(position1, 1, meta); % gain-corrected voltage trace for position1
    position2 = GainCorrectNI(position2, 1, meta); % gain-corrected voltage trace for position2
    laser = GainCorrectNI(laser, 1, meta);
    pacer = GainCorrectNI(pacer, 1, meta);
    
    clearvars temp*
    save('gainCorrectRawTraces', 'camTrig', 'trStart', 'encodeA', 'encodeB', 'reward', 'trEnd', 'lick', 'position1', 'position2', 'laser', 'pacer')
end

%% %% spot the trial-by-trial and all-trials behavioral csv files
% 
% behFilePath = dir(fullfile(p.Results.filePath,'20*')); % dir where the trial-by-trial behavioral csv files are saved
% allTrialCsv = dir(fullfile(behFilePath.folder,behFilePath.name,'trials.csv')); % all trial file
% if length(allTrialCsv)==1
%     trialsFileName = fullfile(allTrialCsv.folder,allTrialCsv.name);
%     trialsCsv = readtable(trialsFileName);
% else
%     error('More than one trials.csv file detected!')
% end

%% task event detection
if ~isempty(dir(fullfile(p.Results.filePath,'evtIndices.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'evtIndices.mat'),'trStartIdx','trEndIdx','rwdIdx','evtIdx8k','evtIdx1k') % if there are gaincorrectedrawtraces already saved, just load them
else
    % detect trialStart, trialEnd
    [trStartIdx,~,~] = detecteventbythreshold(trStart, 8000, 3000, 'stdFactor', 5, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
    fprintf('completed trial start detection!\n');
    [trEndIdx,~,~]   = detecteventbythreshold(trEnd, 8000, 3000, 'stdFactor', 5, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial End
    fprintf('completed trial end detection!\n');
    
    if length(trStartIdx)==length(trEndIdx)
        if sum(trEndIdx - trStartIdx < 0)>=1
            error('Trial End and Start indices do not make sense!')
        end
    elseif length(trStartIdx)-length(trEndIdx)==1
           [trStartIdx] = sort_idx1_by_idx2(trStartIdx, trEndIdx);  
    elseif length(trEndIdx)-length(trStartIdx)==1
           [trEndIdx] = sort_idx2_by_idx1(trStartIdx, trEndIdx); 
    else
        error('Trial End and Start indices do not make sense!')
    end
    
    %% keep detecting events from the nidq file
    % detect reward deliveries
    rwdIdx     = detecteventbythreshold(reward, 8000, 50, 'stdFactor',1, 'plotRez',false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse',true);  % reward
    fprintf('completed reward detection!\n');
    
    % detect camTriggers
    [camTrigRiseIdx, camTrigFallIdx, camPulseTrainIdx] = detecteventbythreshold(camTrig, 8000, 1, 'stdFactor', 1, 'plotRez',false, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', true); % camera trigger
    %camTrigRiseIdx1ms = round(camTrigRiseIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
    %camTrigFallIdx1ms = round(camTrigFallIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
    fprintf('completed camera trigger detection!\n');
    
    % detect laser
    if p.Results.laserUsed
        % detect all laser (laser delivered during the task)
        [evtIdx8k.laserRiseIdx, evtIdx8k.laserFallIdx, evtIdx8k.lsaerPulseTrainIdx] = detecteventbythreshold(laser, 8000, 10, 'stdFactor',2, 'plotRez',true, 'chunkPulses', true, 'chunkInterval', 5000, 'correctLongPulse', true);
    end
    
    % detect position1
    position1Idx     = detecteventbythreshold(position1, 8000, 50, 'stdFactor',1, 'plotRez', true, 'chunkPulses', false, 'detectLater', 0, 'correctLongPulse',true);  % reward
    fprintf('completed position1 detection!\n');
    
    % detect position2
    position2Idx     = detecteventbythreshold(position2, 8000, 50, 'stdFactor',1, 'plotRez', true, 'chunkPulses', false, 'detectLater', 0, 'correctLongPulse',true);  % reward
    fprintf('completed position2 detection!\n');
    
    % detect position3
    position3Idx     = detecteventbythreshold(position3, 8000, 50, 'stdFactor',2, 'plotRez', true, 'chunkPulses', false, 'detectLater', 0, 'correctLongPulse',true);  % reward
    fprintf('completed position3 detection!\n');
    
    % detect pacer
    position3Idx     = detecteventbythreshold(position3, 8000, 50, 'stdFactor',2, 'plotRez', true, 'chunkPulses', false, 'detectLater', 0, 'correctLongPulse',true);  % reward
    fprintf('completed position3 detection!\n');
    
    % store evt indices
    evtIdx8k.trStartIdx = trStartIdx;
    evtIdx8k.trEndIdx = trEndIdx;
    evtIdx8k.rwdIdx  = rwdIdx+nSamp*(p.Results.rewardDelay/1000); % correct for the delay
    evtIdx8k.camTrigRiseIdx = camTrigRiseIdx;
    evtIdx8k.camTrigFallIdx = camTrigFallIdx;
    evtIdx8k.camPulseTrainIdx = camPulseTrainIdx;
    
    evtIdx1k.trStartIdx = round(trStartIdx./8);
    evtIdx1k.trEndIdx = round(trEndIdx./8);
    evtIdx1k.rwdIdx  = round(rwdIdx./8)+p.Results.rewardDelay;
    evtIdx1k.camTrigRiseIdx = round(camTrigRiseIdx./8);
    evtIdx1k.camTrigFallIdx = round(camTrigFallIdx./8);
    evtIdx1k.camPulseTrainIdx = camPulseTrainIdx; % pulse Train Id
    
    save('evtIndices','trStartIdx','trEndIdx','rwdIdx','lickIdx','evtIdx8k','evtIdx1k')
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

jsTime8k = struct;
jsTime1k = struct;

% get the trial information reconstructed based off of trialsCsv, trStartIdx, and rwdIdx
[trialInfo] = reconstructAssay( trialsCsv, trStartIdx, rwdIdx );

for t = 1:length(trStartIdx) % increment trials
    if ~isempty(find(trEndIdx>trStartIdx(t),1)) % if there's a trEnd
        jsTime8k(t).trStart = trStartIdx(t); % trStart detected by the go-cue onset
        % redefine the trial start as the joystick in position timepoint (trJsReadyTime) by examining the baseline encoder data,
        % as the trStart defined with the cue tone onset doesn't align perfectly with the joystick being in the start position.
        jsTime8k(t).trEnd = trEndIdx(find(trEndIdx > trStartIdx(t),1,'first'));  % get the trEnd time point
        
        trBaseRange = trStartIdx(t)-2*nSamp:trStartIdx(t)+(nSamp/10/2)-1; % 2 sec baseline, just take 2 sec before the cue onset, with padding on the righthand side corresponding to 50 ms
        trBasePinA = encodeA(trBaseRange); % pin A this trial baseline
        trBasePinB = encodeB(trBaseRange); % pin B this trial baseline
        biTrBasePinA = double(trBasePinA>pinAmidPoint); % binarize the pinA signal
        biTrBasePinB = double(trBasePinB>pinBmidPoint); % binarize the pinB signal
        
        [trBaseJsState, ~] = readStepperEncoder(biTrBasePinA, biTrBasePinB); % read out the Js position from the binarized quadrature encoded signals
        [trBaseLatestStillTime, jsTime1k(t).baseJsTrajmm, jsTime1k(t).baseSmJsVel, jsTime1k(t).basePeriodicAbsVelSum] = findJsReadyPt2(cumsum(trBaseJsState), nSamp); % find the Js ready time by detecting a still point in the later portion of the baseline period
        
        jsTime8k(t).pull_threshold = trialInfo.pullThreshold(t); % pull threshold
        jsTime8k(t).pull_torque = trialInfo.pullTq(t); % pull torque
        jsTime8k(t).reachP1 = trialInfo.reachPos1(t); % reach position 1
        
        if isnan(trBaseLatestStillTime)
            jsTime8k(t).trJsReady = nan;
            jsTime8k(t).trJsTraj = nan;
            jsTime8k(t).dctrJsTraj = nan;
            jsTime8k(t).rewardT = nan;
            % determine if this trial led to reward
            if ~isempty(find(abs(rwdIdx-jsTime8k(t).trEnd)<=nSamp,1)) % in case there's reward delivery within 1-sec window relative to the trial end
                jsTime8k(t).rewarded = true;
            else
                jsTime8k(t).rewarded = false;
            end
            
            % determine if a laser stim was delivered in this trial
            if p.Results.laserUsed % assign stimLasers to corresponding trials
                if t == 1
                    tempStim = find(evtIdx8k.stimLaserRiseIdx<jsTime8k(t).trEnd);
                    if ~isempty(tempStim)
                        jsTime8k(t).stimLaserOn  = evtIdx8k.stimLaserRiseIdx(tempStim);
                        jsTime8k(t).stimLaserOff = evtIdx8k.stimLaserFallIdx(tempStim);
                    else
                        jsTime8k(t).stimLaserOn  = NaN;
                        jsTime8k(t).stimLaserOff = NaN;
                    end
                    
                elseif t > 1
                    tempStim = find(evtIdx8k.stimLaserRiseIdx<jsTime8k(t).trEnd & evtIdx8k.stimLaserRiseIdx>jsTime8k(t-1).trEnd,1);
                    if ~isempty(tempStim)
                        jsTime8k(t).stimLaserOn  = evtIdx8k.stimLaserRiseIdx(tempStim);
                        jsTime8k(t).stimLaserOff = evtIdx8k.stimLaserFallIdx(tempStim);
                    else
                        jsTime8k(t).stimLaserOn  = NaN;
                        jsTime8k(t).stimLaserOff = NaN;
                    end
                end
            end
            
            if p.Results.plaserUsed
                if t == 1
                    tempPlaser = find(evtIdx8k.plaserRiseIdx<jsTime8k(t).trEnd);
                    if ~isempty(tempPlaser)
                        jsTime8k(t).pLaserOn  = evtIdx8k.plaserRiseIdx(tempPlaser);
                        jsTime8k(t).pLaserOff = evtIdx8k.plaserFallIdx(tempPlaser);
                    else
                        jsTime8k(t).pLaserOn  = NaN;
                        jsTime8k(t).pLaserOff = NaN;
                    end
                    
                elseif t > 1
                    tempPlaser = find(evtIdx8k.plaserRiseIdx<jsTime8k(t).trEnd & evtIdx8k.plaserRiseIdx>jsTime8k(t-1).trEnd,1);
                    if ~isempty(tempPlaser)
                        jsTime8k(t).pLaserOn  = evtIdx8k.plaserRiseIdx(tempPlaser);
                        jsTime8k(t).pLaserOff = evtIdx8k.plaserFallIdx(tempPlaser);
                    else
                        jsTime8k(t).pLaserOn  = NaN;
                        jsTime8k(t).pLaserOff = NaN;
                    end
                end
            end
            
            jsTime8k(t).trialType = nan; % trialType
            jsTime1k(t).movKins = nan; % movement kinematics assign nan
            
        elseif trStartIdx(t)-2*nSamp+trBaseLatestStillTime <= jsTime8k(t).trStart+(nSamp/10/2) % if there's a joystick still point within 50 ms from the trStart signal
            jsTime8k(t).trJsReady = trStartIdx(t)-2*nSamp+trBaseLatestStillTime; % get the trJsReady time
            
            trRange = jsTime8k(t).trJsReady:jsTime8k(t).trEnd; % trial range aligned to the Js ready time
            
            trPinA = encodeA(trRange); % pin A this trial aligned to the Js ready
            trPinB = encodeB(trRange); % pin A this trial aligned to the Js ready
            
            biTrPinA = double(trPinA>pinAmidPoint); % binarize the encoder signal
            biTrPinB = double(trPinB>pinBmidPoint); % binarize the encoder signal
            [trJsState, ~] = readStepperEncoder(biTrPinA, biTrPinB); % read out the Js position from the binarized quadrature encoded signals
            jsTime8k(t).trJsTraj = cumsum(trJsState);  % get the cumulative joystick trajectory for this trial
            jsTime8k(t).dctrJsTraj = decimate(jsTime8k(t).trJsTraj ,round(nSamp/1000)); % decimate the Traj (downsampling the 8kHz data at 1kHz, so the time resoluation to be 1ms)
            
            % to use sg filter get sgfiltFramelen
            if length(jsTime8k(t).dctrJsTraj) >= p.Results.sgfiltFramelen
                sgfiltFramelen = p.Results.sgfiltFramelen;
            elseif length(jsTime8k(t).dctrJsTraj) < p.Results.sgfiltFramelen
                if mod(length(jsTime8k(t).dctrJsTraj),2)==0
                    sgfiltFramelen = length(jsTime8k(t).dctrJsTraj)-1; % the frame length for sg filter needs to be an odd number
                else
                    sgfiltFramelen = length(jsTime8k(t).dctrJsTraj);
                end
            end
            
            tempSmdctrJsTraj = sgolayfilt(jsTime8k(t).dctrJsTraj,3,sgfiltFramelen);
            %tempSmdctrJsTrajmm = tempSmdctrJsTraj*jsDistCoeff; % convert the Js traj into mm
            
            % determine if the trial got rewarded
            if ~isempty(find(abs(rwdIdx-jsTime8k(t).trEnd)<=nSamp,1)) % in case there's reward delivery within 1-sec window relative to the trial end
                jsTime8k(t).rewarded = true;
            else
                jsTime8k(t).rewarded = false;
            end
            
            % determine if a laser stim was delivered in this trial
            if p.Results.laserUsed % assign stimLasers to corresponding trials
                if t == 1
                    tempStim = find(evtIdx8k.stimLaserRiseIdx<jsTime8k(t).trEnd);
                    if ~isempty(tempStim)
                        jsTime8k(t).stimLaserOn  = evtIdx8k.stimLaserRiseIdx(tempStim);
                        jsTime8k(t).stimLaserOff = evtIdx8k.stimLaserFallIdx(tempStim);
                    else
                        jsTime8k(t).stimLaserOn  = NaN;
                        jsTime8k(t).stimLaserOff = NaN;
                    end
                elseif t > 1
                    tempStim = find(evtIdx8k.stimLaserRiseIdx<jsTime8k(t).trEnd & evtIdx8k.stimLaserRiseIdx>jsTime8k(t-1).trEnd,1);
                    if ~isempty(tempStim)
                        jsTime8k(t).stimLaserOn  = evtIdx8k.stimLaserRiseIdx(tempStim);
                        jsTime8k(t).stimLaserOff = evtIdx8k.stimLaserFallIdx(tempStim);
                    else
                        jsTime8k(t).stimLaserOn  = NaN;
                        jsTime8k(t).stimLaserOff = NaN;
                    end
                end
            end
            
            if p.Results.plaserUsed
                if t == 1
                    tempPlaser = find(evtIdx8k.plaserRiseIdx<jsTime8k(t).trEnd);
                    if ~isempty(tempPlaser)
                        jsTime8k(t).pLaserOn  = evtIdx8k.plaserRiseIdx(tempPlaser);
                        jsTime8k(t).pLaserOff = evtIdx8k.plaserFallIdx(tempPlaser);
                    else
                        jsTime8k(t).pLaserOn  = NaN;
                        jsTime8k(t).pLaserOff = NaN;
                    end
                    
                elseif t > 1
                    tempPlaser = find(evtIdx8k.plaserRiseIdx<jsTime8k(t).trEnd & evtIdx8k.plaserRiseIdx>jsTime8k(t-1).trEnd,1);
                    if ~isempty(tempPlaser)
                        jsTime8k(t).pLaserOn  = evtIdx8k.plaserRiseIdx(tempPlaser);
                        jsTime8k(t).pLaserOff = evtIdx8k.plaserFallIdx(tempPlaser);
                    else
                        jsTime8k(t).pLaserOn  = NaN;
                        jsTime8k(t).pLaserOff = NaN;
                    end
                end
            end
            
            % classify the trial ('sp': successfull pull, 'ps': push, 'pm': premature pull, 'to': timeout, 'nn': not identified)
            if jsTime8k(t).rewarded % if rewarded
                if ~isempty(find(tempSmdctrJsTraj<jsTime8k(t).pull_threshold,1)) % check the negative threshold crossing
                    jsTime8k(t).trialType = 'sp'; % successful pull
                    jsTime8k(t).rewardT = jsTime8k(t).trEnd+nSamp*(p.Results.rewardDelay/1000); % reward if delivered must be 1000ms (or some other delay) after the trial end
                else % this should be an premature pull which has been accidentally rewarded
                    jsTime8k(t).trialType = nan; % just throw out these trials for now
                    jsTime8k(t).rewardT = nan;
                end
            else % if not rewarded
                jsTime8k(t).rewardT = nan;
                if length(jsTime8k(t).dctrJsTraj)>p.Results.trialTimeout-100 % timeout
                    jsTime8k(t).trialType = 'to';
                elseif isempty(find(jsTime8k(t).dctrJsTraj<jsTime8k(t).pull_threshold,1)) && ~isempty(find(jsTime8k(t).dctrJsTraj>50,1)) % push (no pull beyond the pull threshold && push beyond a certain threshold)
                    jsTime8k(t).trialType = 'ps';
                elseif ~isempty(find(jsTime8k(t).dctrJsTraj<jsTime8k(t).pull_threshold,1)) % pull (unrewarded)
                    impullThresCross = find(jsTime8k(t).dctrJsTraj<jsTime8k(t).pull_threshold,1); % premature pull threshold crossing point
                    if ~isempty(find(jsTime8k(t).dctrJsTraj(impullThresCross:end)>0,1))
                        jsTime8k(t).trialType = 'pmpp'; % premature pull and push
                    else
                        jsTime8k(t).trialType = 'pm'; % premature pull
                    end
                else
                    jsTime8k(t).trialType = 'nn';
                end
            end
            
            % trial-by-trial reach kinematics analysis
            tempMass = p.Results.meanMass(2,jsTime8k(t).pull_torque==p.Results.meanMass(1,:));
            [ jsTime1k(t).movKins ] = jsReachKinematics( jsTime8k(t).dctrJsTraj, jsTime8k(t).pull_threshold, jsTime8k(t).trialType, tempMass, sgfiltFramelen );
            %figure; plot(interp1(1:length(jsTime8k(t).csvTraj),jsTime8k(t).csvTraj,1:length(jsTime8k(t).dcsmtrJsTraj)));
        end
    else % if there's no trialEnd left
    end
    fprintf('completed trial #%d\n', t);
end
clearvars temp* t

% transfer relevant information to jsTime1k struct from jsTime8k
n2cTrStart = num2cell(round([jsTime8k(:).trStart]./8)); [jsTime1k.trStart] = n2cTrStart{:};
n2cTrEnd   = num2cell(round([jsTime8k(:).trEnd]./8));   [jsTime1k.trEnd]   = n2cTrEnd{:};
n2cTrJsReady = num2cell(round([jsTime8k(:).trJsReady]./8)); [jsTime1k.trJsReady] = n2cTrJsReady{:};
n2cRewarded = num2cell([jsTime8k(:).rewarded]); [jsTime1k.rewarded] = n2cRewarded{:};
n2cPull_threshold = num2cell([jsTime8k(:).pull_threshold]); [jsTime1k.pull_threshold] = n2cPull_threshold{:};
n2cPull_torque = num2cell([jsTime8k(:).pull_torque]); [jsTime1k.pull_torque] = n2cPull_torque{:};
n2cReachP1 = num2cell([jsTime8k(:).reachP1]); [jsTime1k.reachP1] = n2cReachP1{:};
n2cRewardT = num2cell(round([jsTime8k(:).rewardT]./8)); [jsTime1k.rewardT] = n2cRewardT{:};
if p.Results.laserUsed
    n2cStimLaserOn = num2cell(round([jsTime8k(:).stimLaserOn]./8)); [jsTime1k.stimLaserOn] = n2cStimLaserOn{:};
    n2cStimLaserOff = num2cell(round([jsTime8k(:).stimLaserOff]./8)); [jsTime1k.stimLaserOff] = n2cStimLaserOff{:};
end
if p.Results.plaserUsed
    n2cPLaserOn = num2cell(round([jsTime8k(:).pLaserOn]./8)); [jsTime1k.pLaserOn] = n2cPLaserOn{:};
    n2cPLaserOff = num2cell(round([jsTime8k(:).pLaserOff]./8)); [jsTime1k.pLaserOff] = n2cPLaserOff{:};
end

if exist('imecFiles','var')==1
    if ~isempty(imecFiles)
        if length(imecFiles)==1
            if isequal(length(trStartIdx),length(trStartImec)) && isequal(length(trEndIdx),length(trEndImec))
                n2cTrStartImec = num2cell(trStartImec); [jsTime1k.trStartImec] = n2cTrStartImec{:};
                n2cTrEndImec = num2cell(trEndImec); [jsTime1k.trEndImec] = n2cTrEndImec{:};
            end
        end
    end
end

[jsTime1k.trialType] = jsTime8k(:).trialType;

% generate a plot to inspect a certain trial
%movKinsPlot(jsTime1k(t).movKins);

%% Save relevant BehVariables
cd(p.Results.filePath)
save('BehVariablesSwitch', 'jsTime1k', 'jsTime8k', 'evtIdx8k', 'evtIdx1k', 'p', 'trialsCsv', 'trialInfo', 'tbytCsvList') % append the position/velocity data variables

end

function [trStartImec, trEndImec] = saveImecEvent(data_file, event_data, meta)
clearvars trStart trEnd

% 0: trial start (note that matlab is 1-based)
% 7: trial end

trStartEvt = double(bitget(event_data, 1, 'uint16')); %[0; diff(double(bitget(event_data, 1, 'uint16')))];
[trStartIdx,~,~] = detecteventbythreshold(trStartEvt', meta.imSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
trStartImec = trStartIdx./meta.imSampRate*1000; % time in msec

trEndEvt = double(bitget(event_data, 8, 'uint16'));   %[0; diff(double(bitget(event_data, 1, 'uint16')))];
[trEndIdx,~,~] = detecteventbythreshold(trEndEvt', meta.imSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial Start
trEndImec = trEndIdx./meta.imSampRate*1000; % time in msec

disp(['Saving IMEC data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'trStartImec', 'trEndImec', '-append');
else
    save(data_file, 'trStartImec', 'trEndImec');
end
end

function [correct_idx2] = sort_idx2_by_idx1(idx1, idx2)
        % sort index 2 relative to index 1
        correct_idx2 = nan(1, length(idx1)); 
        for j = 1:length(idx1) 
            correct_idx = find(idx1(j)<idx2, 1, 'first'); 
            correct_idx2(1,j) = idx2(correct_idx); 
        end
end

function [correct_idx1] = sort_idx1_by_idx2(idx1, idx2)
        % sort index 1 relative to index 2
        correct_idx1 = nan(1, length(idx1)); 
        for j = 1:length(idx2) 
            correct_idx = find(idx1<idx2(j), 1, 'last'); 
            correct_idx1(1,j) = idx1(correct_idx); 
        end
end








