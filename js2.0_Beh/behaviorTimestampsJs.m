function behaviorTimestampsJs(p)
%behaviorTimestamps

%addpath(genpath(''))
%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718';
%p = parse_input_Js(filePath, varargin ); % parse input
%p = parse_input_Js(filePath, {} ); % use this line instead when running line-by-line
% To convert datenum to datetime use, datetime(vFronFileStartDatenum(tempVFronI),'ConvertFrom','datenum')

cd(p.Results.filePath)

% check the folder whether there's rez file already
if ~isempty(dir(fullfile(p.Results.filePath,'BehVariablesJs.mat'))) % if the BehVariablesJs.mat file already exists in the filePath
    answer = questdlg('BehVariablesJs.mat already exists, Would you like to replace it?','Choice','Replace','Cancel','Cancel');
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
nSamp = round(SampRate(meta),0); % sampling rate (default: 25kHz)
%totalTimeSecsMeta = str2double(meta.fileTimeSecs); % total duration of file in seconds

if ~isempty(dir(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'), 'lick', 'trStart', 'reward', 'trEnd', 'camTrig', 'encodeA', 'encodeB', 'laser', 'plaser' ) % if there are gaincorrectedrawtraces already saved, just load them
else
    % Specify the relevant behavioral channel numbers
    trStartCh  = channels{1}+p.Results.trStartCh; % ch# for trial start
    camTrigCh  = channels{1}+p.Results.camTrigCh; % ch# for camera trigger
    rewardCh   = channels{1}+p.Results.rewardCh;  % ch# for reward delivery
    trEndCh    = channels{1}+p.Results.trEndCh;   % ch# for trial end (either by successful pull or error/timeout)
    encodeACh  = channels{1}+p.Results.encodeACh; % ch# for stepper direction
    encodeBCh = channels{1}+p.Results.encodeBCh; % ch# for stepper steps
    laserCh   = channels{1}+p.Results.laserCh;   % ch# for laser triggers
    plaserCh  = channels{1}+p.Results.plaserCh;  % ch# for pseudolaser triggers
    lickCh = channels{1}+p.Results.lickCh;      % ch# for lick detect
    
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
    trStart   = zeros(1,floor(totalTimeSecs*25000));
    camTrig   = zeros(1,floor(totalTimeSecs*25000));
    reward    = zeros(1,floor(totalTimeSecs*25000));
    trEnd     = zeros(1,floor(totalTimeSecs*25000));
    encodeA   = zeros(1,floor(totalTimeSecs*25000));
    encodeB   = zeros(1,floor(totalTimeSecs*25000));
    laser     = zeros(1,floor(totalTimeSecs*25000));
    plaser    = zeros(1,floor(totalTimeSecs*25000));
    lick      = zeros(1,floor(totalTimeSecs*25000));
    
    for k = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        tempDataArray = ReadBin(k*nSamp, nSamp, meta, binName, binFile.folder); % read bin data for each second
        tempTrStart  = tempDataArray(trStartCh,:); % decimate the data
        tempCamTrig  = tempDataArray(camTrigCh,:); % do not decimate for higher temporal resolution
        tempReward   = tempDataArray(rewardCh,:);
        tempTrEnd    = tempDataArray(trEndCh,:);
        tempEncodeA  = tempDataArray(encodeACh,:); % do not decimate for higher temporal resolution
        tempEncodeB  = tempDataArray(encodeBCh,:); % do not decimate for higher temporal resolution
        tempLick     = tempDataArray(lickCh,:);
        tempLaser    = tempDataArray(laserCh,:);
        tempPLaser   = tempDataArray(plaserCh,:);
        %tempSync = tempDataArray(3,:); 
        
        trStart(1,k*25000+1:(k+1)*25000) = tempTrStart; % accumulated the decimated data second-by-second
        camTrig(1,k*25000+1:(k+1)*25000) = tempCamTrig;
        reward(1,k*25000+1:(k+1)*25000) = tempReward;
        trEnd(1,k*25000+1:(k+1)*25000)  = tempTrEnd;
        encodeA(1,k*25000+1:(k+1)*25000) = tempEncodeA;
        encodeB(1,k*25000+1:(k+1)*25000)  = tempEncodeB;
        lick(1,k*25000+1:(k+1)*25000) = tempLick;
        laser(1,k*25000+1:(k+1)*25000) = tempLaser;
        plaser(1,k*25000+1:(k+1)*25000) = tempPLaser;
        %sync(1,k*25000+1:(k+1)*25000) = tempSync; 
        
        fprintf('processed %d\n', k+1)
    end
    clearvars k
    
    % Gain correction for channnels of interest
    if strcmp(meta.typeThis, 'imec') % in case recording via imec
        trStart = GainCorrectIM(trStart, 1, meta); % gain-corrected voltage trace for trStart
        camTrig = GainCorrectIM(camTrig, 1, meta); % gain-corrected voltage trace for camTrig
        reward = GainCorrectIM(reward, 1, meta);   % gain-corrected voltage trace for reward
        trEnd  = GainCorrectIM(trEnd, 1, meta);    % gain-corrected voltage trace for trEnd
        encodeA = GainCorrectIM(encodeA, 1, meta); % gain-corrected voltage trace for encodeA
        encodeB = GainCorrectIM(encodeB, 1, meta); % gain-corrected voltage trace for encodeB
        lick = GainCorrectIM(lick, 1, meta); % gain-corrected voltage trace for lick
        laser = GainCorrectIM(laser, 1, meta);
        plaser = GainCorrectIM(plaser, 1, meta);
    else    % in case of recording via NI board
        trStart = GainCorrectNI(trStart, 1, meta); % gain-corrected voltage trace for trStart
        camTrig = GainCorrectNI(camTrig, 1, meta); % gain-corrected voltage trace for camTrig
        reward = GainCorrectNI(reward, 1, meta);   % gain-corrected voltage trace for reward
        trEnd  = GainCorrectNI(trEnd, 1, meta);    % gain-corrected voltage trace for trEnd
        encodeA = GainCorrectNI(encodeA, 1, meta); % gain-corrected voltage trace for encodeA
        encodeB = GainCorrectNI(encodeB, 1, meta); % gain-corrected voltage trace for encodeB
        lick = GainCorrectNI(lick, 1, meta); % gain-corrected voltage trace for lick
        laser = GainCorrectNI(laser, 1, meta);
        plaser = GainCorrectNI(plaser, 1, meta);
    end
    clearvars temp*
    save('gainCorrectRawTraces', 'trStart', 'camTrig', 'reward', 'trEnd', 'encodeA', 'encodeB', 'lick', 'laser', 'plaser')
end

%% spot the trial-by-trial and all-trials behavioral csv files
if isempty(dir(fullfile(p.Results.filePath,'20*')))
    error('Cannot find the trial-by-trial behavior data csv files!')
end

behFilePath = dir(fullfile(p.Results.filePath,'20*')); % dir where the trial-by-trial behavioral csv files are saved
tbytCsvList = dir(fullfile(behFilePath.folder,behFilePath.name,'trial_*'));    % trial-by-trial files
allTrialCsv = dir(fullfile(behFilePath.folder,behFilePath.name,'trials.csv')); % all trial file
if length(allTrialCsv)==1
    trialsFileName = fullfile(allTrialCsv.folder,allTrialCsv.name);
    trialsCsv = readtable(trialsFileName);
else
    error('More than one trials.csv file detected!')
end

[~,tbytCsvdateSort] = sort(datenum({tbytCsvList(:).date}, 'dd-mmm-yyyy hh:MM:ss'), 1, 'ascend'); % sorted fileList

%% task event detection
if ~isempty(dir(fullfile(p.Results.filePath,'evtIndices.mat'))) && p.Results.reReadBin==false % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'evtIndices.mat'),'trStartIdx','trEndIdx','rwdIdx','lickIdx','evtIdx25k','evtIdx1k') % if there are gaincorrectedrawtraces already saved, just load them
else
    [trStartIdx,~,~] = detecteventbythreshold(trStart, 25000, 3000, 'stdFactor', 5, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
    fprintf('completed trial start detection!');
    [trEndIdx,~,~]   = detecteventbythreshold(trEnd, 25000, 3000, 'stdFactor', 5, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial End
    fprintf('completed trial end detection!');
    
    if length(trStartIdx)==length(trEndIdx)
        if ~unique(trEndIdx - trStartIdx>0)
            error('Trial End and Start indices do not make sense!')
        end
    elseif length(trStartIdx)-length(trEndIdx)==1
        if ~unique(trEndIdx - trStartIdx(1,1:length(trEndIdx))>0)
            error('Trial End and Start indices do not make sense!')
        end
    elseif length(trEndIdx)-length(trStartIdx)==1
        if ~unique(trEndIdx(1,1:length(trStartIdx)) - trStartIdx>0)
            error('Trial End and Start indices do not make sense!')
        end
    elseif length(trStartIdx)==size(trialsCsv,1)
        newTrEndIdx = sortTrStartTrEnd(trStartIdx, trEndIdx, trialsCsv, p.Results.trialTimeout);
        trEndIdx = newTrEndIdx; 
    else
        error('Trial End and Start indices do not make sense!')
    end
    
    %% see if there's an imec data associated with the current file
    if contains(p.Results.filePath,'ni','IgnoreCase',true)
        subFilePath = strsplit(p.Results.filePath,'\');
        imecSearchStartPath = fullfile(subFilePath{1,1:end-1});
    else
        imecSearchStartPath = p.Results.filePath;
    end
    subPath = strsplit(genpath(imecSearchStartPath),';');
    imecFiles = [];
    for f = 1:length(subPath)-1
        imecFiles = [imecFiles; dir(fullfile(subPath{f},'*.imec.ap.bin'))];
    end
    clearvars f
    
    if length(imecFiles)==1
        % data file name
        disp(['======== ', imecFiles(1).name, ' ========']);
        % read the imec bin file and detect trStart and trEnd events
        [event_imec, meta_imec] = readEventBin(fullfile(imecFiles(1).folder,imecFiles(1).name)); %save.readEventBin(imec_file{i_bin});
        imecEvtFile = fullfile(imecFiles(1).folder,'imecEvtData.mat');
        [trStartImec, trEndImec] = saveImecEvent(imecEvtFile, event_imec, meta_imec);
        clear event_imec
        
    elseif length(imecFiles)>1
        error('Multiple imec.ap.bin files were detected!')
    end
    
    %% keep detecting events from the nidq file
    % detect reward deliveries
    rwdIdx     = detecteventbythreshold(reward, 25000, 50, 'stdFactor',1, 'plotRez',false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse',true);  % reward
    fprintf('completed reward detection!');
    
    % detect licks
    lickIdx = detecteventbythreshold(decimate(lick,nSamp/1000), 1000, 20, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false); % lick, detect licks after downsampling (there seems to be some denoising effect with decimation, see /Volumes/RAID2/parkj/NeuralData/js2.0/WR25/101718/lickSignalExample_1kHzVS25kHz.fig as an example)
    %deciLick = decimate(lick, nSamp/1000); intDeciLick = interp1(1:length(deciLick), deciLick, linespace(1, length(deciLick), length(lick)));
    %plot(lick); hold on; plot(intDeciLick); hold off
    fprintf('completed lick detection!');
    
    % detect camTriggers
    [camTrigRiseIdx, camTrigFallIdx, camPulseTrainIdx] = detecteventbythreshold(camTrig, 25000, 2, 'stdFactor', 1, 'plotRez',false, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', true); % camera trigger
    %camTrigRiseIdx1ms = round(camTrigRiseIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
    %camTrigFallIdx1ms = round(camTrigFallIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
    fprintf('completed camera trigger detection!');
    
    % detect laser
    if p.Results.laserUsed
        % detect all laser (laser delivered during the task)
        [evtIdx25k.laserRiseIdx, evtIdx25k.laserFallIdx] = detecteventbythreshold(laser, 25000, 50, 'stdFactor',2, 'plotRez',false, 'chunkPulses', false);
        evtIdx25k.stimLaserRiseIdx = evtIdx25k.laserRiseIdx(evtIdx25k.laserRiseIdx<trEndIdx(end));
        evtIdx25k.stimLaserFallIdx = evtIdx25k.laserFallIdx(evtIdx25k.laserFallIdx<trEndIdx(end));
        evtIdx1k.stimLaserRiseIdx = round(evtIdx25k.stimLaserRiseIdx./25);
        evtIdx1k.stimLaserFallIdx = round(evtIdx25k.stimLaserFallIdx./25);
        if p.Results.tagLaserUsed
            evtIdx25k.tagLaserRiseIdx = evtIdx25k.laserRiseIdx(end-p.Results.numbTagLaser+1:end);
            evtIdx25k.tagLaserFallIdx = evtIdx25k.laserFallIdx(end-p.Results.numbTagLaser+1:end);
            evtIdx1k.tagLaserRiseIdx = round(evtIdx25k.tagLaserRiseIdx./25);
            evtIdx1k.tagLaserFallIdx = round(evtIdx25k.tagLaserFallIdx./25);
        end
    end
    
    if p.Results.plaserUsed
        [evtIdx25k.plaserRiseIdx, evtIdx25k.plaserFallIdx] = detecteventbythreshold(plaser, 25000, 50, 'stdFactor',1, 'plotRez',false, 'chunkPulses', false);
        evtIdx1k.plaserRiseIdx = round(evtIdx25k.plaserRiseIdx./25);
        evtIdx1k.plaserFallIdx = round(evtIdx25k.plaserFallIdx./25);
    end
    
    % store evt indices
    evtIdx25k.trStartIdx = trStartIdx;
    evtIdx25k.trEndIdx = trEndIdx;
    evtIdx25k.rwdIdx  = rwdIdx+nSamp*(p.Results.rewardDelay/1000); % correct for the delay
    %evtIdx25k.lickIdx = lickIdx;
    evtIdx25k.camTrigRiseIdx = camTrigRiseIdx;
    evtIdx25k.camTrigFallIdx = camTrigFallIdx;
    evtIdx25k.camPulseTrainIdx = camPulseTrainIdx;
    
    evtIdx1k.trStartIdx = round(trStartIdx./25);
    evtIdx1k.trEndIdx = round(trEndIdx./25);
    
    if ~isempty(imecFiles)
        if length(imecFiles)==1
            if isequal(length(trStartIdx),length(trStartImec)) && isequal(length(trEndIdx),length(trEndImec))
                evtIdx1k.trStartImec = trStartImec;
                evtIdx1k.trEndImec = trEndImec;
            else
                warning('The # of trStart or trEnd files detected from nidq and imec differs!')
            end
        end
    end
    
    evtIdx1k.rwdIdx  = round(rwdIdx./25)+p.Results.rewardDelay;
    evtIdx1k.lickIdx = round(lickIdx);
    evtIdx1k.camTrigRiseIdx = round(camTrigRiseIdx./25);
    evtIdx1k.camTrigFallIdx = round(camTrigFallIdx./25);
    evtIdx1k.camPulseTrainIdx = camPulseTrainIdx; % pulse Train Id
    
    save('evtIndices','trStartIdx','trEndIdx','rwdIdx','lickIdx','evtIdx25k','evtIdx1k')
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
[trialInfo] = reconstructAssay( trialsCsv, trStartIdx, rwdIdx );

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
            
            % determine if a laser stim was delivered in this trial
            if p.Results.laserUsed % assign stimLasers to corresponding trials
                if t == 1
                    tempStim = find(evtIdx25k.stimLaserRiseIdx<jsTime25k(t).trEnd);
                    if ~isempty(tempStim)
                        jsTime25k(t).stimLaserOn  = evtIdx25k.stimLaserRiseIdx(tempStim);
                        jsTime25k(t).stimLaserOff = evtIdx25k.stimLaserFallIdx(tempStim);
                    else
                        jsTime25k(t).stimLaserOn  = NaN;
                        jsTime25k(t).stimLaserOff = NaN;
                    end
                    
                elseif t > 1
                    tempStim = find(evtIdx25k.stimLaserRiseIdx<jsTime25k(t).trEnd & evtIdx25k.stimLaserRiseIdx>jsTime25k(t-1).trEnd,1);
                    if ~isempty(tempStim)
                        jsTime25k(t).stimLaserOn  = evtIdx25k.stimLaserRiseIdx(tempStim);
                        jsTime25k(t).stimLaserOff = evtIdx25k.stimLaserFallIdx(tempStim);
                    else
                        jsTime25k(t).stimLaserOn  = NaN;
                        jsTime25k(t).stimLaserOff = NaN;
                    end
                end
            end
            
            if p.Results.plaserUsed
                if t == 1
                    tempPlaser = find(evtIdx25k.plaserRiseIdx<jsTime25k(t).trEnd);
                    if ~isempty(tempPlaser)
                        jsTime25k(t).pLaserOn  = evtIdx25k.plaserRiseIdx(tempPlaser);
                        jsTime25k(t).pLaserOff = evtIdx25k.plaserFallIdx(tempPlaser);
                    else
                        jsTime25k(t).pLaserOn  = NaN;
                        jsTime25k(t).pLaserOff = NaN;
                    end
                    
                elseif t > 1
                    tempPlaser = find(evtIdx25k.plaserRiseIdx<jsTime25k(t).trEnd & evtIdx25k.plaserRiseIdx>jsTime25k(t-1).trEnd,1);
                    if ~isempty(tempPlaser)
                        jsTime25k(t).pLaserOn  = evtIdx25k.plaserRiseIdx(tempPlaser);
                        jsTime25k(t).pLaserOff = evtIdx25k.plaserFallIdx(tempPlaser);
                    else
                        jsTime25k(t).pLaserOn  = NaN;
                        jsTime25k(t).pLaserOff = NaN;
                    end
                end
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
            
            % determine if a laser stim was delivered in this trial
            if p.Results.laserUsed % assign stimLasers to corresponding trials
                if t == 1
                    tempStim = find(evtIdx25k.stimLaserRiseIdx<jsTime25k(t).trEnd);
                    if ~isempty(tempStim)
                        jsTime25k(t).stimLaserOn  = evtIdx25k.stimLaserRiseIdx(tempStim);
                        jsTime25k(t).stimLaserOff = evtIdx25k.stimLaserFallIdx(tempStim);
                    else
                        jsTime25k(t).stimLaserOn  = NaN;
                        jsTime25k(t).stimLaserOff = NaN;
                    end
                elseif t > 1
                    tempStim = find(evtIdx25k.stimLaserRiseIdx<jsTime25k(t).trEnd & evtIdx25k.stimLaserRiseIdx>jsTime25k(t-1).trEnd,1);
                    if ~isempty(tempStim)
                        jsTime25k(t).stimLaserOn  = evtIdx25k.stimLaserRiseIdx(tempStim);
                        jsTime25k(t).stimLaserOff = evtIdx25k.stimLaserFallIdx(tempStim);
                    else
                        jsTime25k(t).stimLaserOn  = NaN;
                        jsTime25k(t).stimLaserOff = NaN;
                    end
                end
            end
            
            if p.Results.plaserUsed
                if t == 1
                    tempPlaser = find(evtIdx25k.plaserRiseIdx<jsTime25k(t).trEnd);
                    if ~isempty(tempPlaser)
                        jsTime25k(t).pLaserOn  = evtIdx25k.plaserRiseIdx(tempPlaser);
                        jsTime25k(t).pLaserOff = evtIdx25k.plaserFallIdx(tempPlaser);
                    else
                        jsTime25k(t).pLaserOn  = NaN;
                        jsTime25k(t).pLaserOff = NaN;
                    end
                    
                elseif t > 1
                    tempPlaser = find(evtIdx25k.plaserRiseIdx<jsTime25k(t).trEnd & evtIdx25k.plaserRiseIdx>jsTime25k(t-1).trEnd,1);
                    if ~isempty(tempPlaser)
                        jsTime25k(t).pLaserOn  = evtIdx25k.plaserRiseIdx(tempPlaser);
                        jsTime25k(t).pLaserOff = evtIdx25k.plaserFallIdx(tempPlaser);
                    else
                        jsTime25k(t).pLaserOn  = NaN;
                        jsTime25k(t).pLaserOff = NaN;
                    end
                end
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
            tempMass = p.Results.meanMass(2,jsTime25k(t).pull_torque==p.Results.meanMass(1,:));
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
if p.Results.plaserUsed
    n2cPLaserOn = num2cell(round([jsTime25k(:).pLaserOn]./25)); [jsTime1k.pLaserOn] = n2cPLaserOn{:};
    n2cPLaserOff = num2cell(round([jsTime25k(:).pLaserOff]./25)); [jsTime1k.pLaserOff] = n2cPLaserOff{:};
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

[jsTime1k.trialType] = jsTime25k(:).trialType;

% generate a plot to inspect a certain trial
%movKinsPlot(jsTime1k(t).movKins);

%% Save relevant BehVariables
cd(p.Results.filePath)
save('BehVariablesJs', 'jsTime1k', 'jsTime25k', 'evtIdx25k', 'evtIdx1k', 'p', 'trialsCsv', 'trialInfo', 'tbytCsvList') % append the position/velocity data variables

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
