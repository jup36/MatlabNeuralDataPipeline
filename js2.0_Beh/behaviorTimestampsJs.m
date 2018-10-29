function behaviorTimestampsJs(filePath, varargin)
%behaviorTimestamps

%addpath(genpath(''))
%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR28';
p = parse_input_Js(filePath, varargin ); % parse input
%p = parse_input_Js(filePath, {} ); % use this line instead when running line-by-line

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
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Parse the corresponding metafile
meta  = ReadMeta(binName, p.Results.filePath); % get the meta data (structure)
nSamp = SampRate(meta);          % sampling rate (default: 25kHz)
totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds

if ~isempty(dir(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'))) % if the gainCorrectRawTraces.mat file already exists in the filePath
    load(fullfile(p.Results.filePath,'gainCorrectRawTraces.mat'), 'lick', 'trStart', 'reward', 'trEnd', 'camTrig', 'encodeA', 'encodeB' ) % if there are gaincorrectedrawtraces already saved, just load them
else
    
    % Specify the relevant behavioral channel numbers
    trStartCh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.trStartCh; % ch# for trial start
    camTrigCh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.camTrigCh; % ch# for camera trigger
    rewardCh   = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.rewardCh;  % ch# for reward delivery
    trEndCh    = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.trEndCh;   % ch# for trial end (either by successful pull or error/timeout)
    encodeACh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.encodeACh;  % ch# for stepper direction
    encodeBCh = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.encodeBCh; % ch# for stepper steps
    lickCh = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.lickCh;   % ch# for lick detect
    
    % preallocate the behavioral data arrays
    trStart   = zeros(1,floor(totalTimeSecs*25000));  % the time resolution will be 1000Hz (1ms) after decimation
    camTrig   = zeros(1,floor(totalTimeSecs*25000)); % do not decimate
    reward    = zeros(1,floor(totalTimeSecs*25000));  %
    trEnd     = zeros(1,floor(totalTimeSecs*25000));  %
    encodeA   = zeros(1,floor(totalTimeSecs*25000)); % do not decimate
    encodeB   = zeros(1,floor(totalTimeSecs*25000)); % do not decimate
    lick      = zeros(1,floor(totalTimeSecs*25000));  %
    
    for i = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        tempDataArray = ReadBin(i*nSamp, nSamp, meta, binName, p.Results.filePath); % read bin data for each second
        tempTrStart  = tempDataArray(trStartCh,:); % decimate the data
        tempCamTrig  = tempDataArray(camTrigCh,:); % do not decimate for higher temporal resolution
        tempReward   = tempDataArray(rewardCh,:);
        tempTrEnd    = tempDataArray(trEndCh,:);
        tempEncodeA  = tempDataArray(encodeACh,:); % do not decimate for higher temporal resolution
        tempEncodeB  = tempDataArray(encodeBCh,:); % do not decimate for higher temporal resolution
        tempLick     = tempDataArray(lickCh,:);
        
        trStart(1,i*25000+1:(i+1)*25000) = tempTrStart; % accumulated the decimated data second-by-second
        camTrig(1,i*25000+1:(i+1)*25000) = tempCamTrig;
        reward(1,i*25000+1:(i+1)*25000) = tempReward;
        trEnd(1,i*25000+1:(i+1)*25000)  = tempTrEnd;
        encodeA(1,i*25000+1:(i+1)*25000) = tempEncodeA;
        encodeB(1,i*25000+1:(i+1)*25000)  = tempEncodeB;
        lick(1,i*25000+1:(i+1)*25000) = tempLick;
        fprintf('processed %d\n', i+1)
    end
    clearvars i
    
    % Gain correction for channnels of interest
    if strcmp(meta.typeThis, 'imec') % in case recording via imec
        trStart = GainCorrectIM(trStart, 1, meta); % gain-corrected voltage trace for trStart
        camTrig = GainCorrectIM(camTrig, 1, meta); % gain-corrected voltage trace for camTrig
        reward = GainCorrectIM(reward, 1, meta);   % gain-corrected voltage trace for reward
        trEnd  = GainCorrectIM(trEnd, 1, meta);    % gain-corrected voltage trace for trEnd
        encodeA = GainCorrectIM(encodeA, 1, meta); % gain-corrected voltage trace for encodeA
        encodeB = GainCorrectIM(encodeB, 1, meta); % gain-corrected voltage trace for encodeB
        lick = GainCorrectIM(lick, 1, meta); % gain-corrected voltage trace for lick
    else    % in case of recording via NI board
        trStart = GainCorrectNI(trStart, 1, meta); % gain-corrected voltage trace for trStart
        camTrig = GainCorrectNI(camTrig, 1, meta); % gain-corrected voltage trace for camTrig
        reward = GainCorrectNI(reward, 1, meta);   % gain-corrected voltage trace for reward
        trEnd  = GainCorrectNI(trEnd, 1, meta);    % gain-corrected voltage trace for trEnd
        encodeA = GainCorrectNI(encodeA, 1, meta); % gain-corrected voltage trace for encodeA
        encodeB = GainCorrectNI(encodeB, 1, meta); % gain-corrected voltage trace for encodeB
        lick = GainCorrectNI(lick, 1, meta); % gain-corrected voltage trace for lick
    end
    clearvars temp*
    save('gainCorrectRawTraces', 'trStart', 'camTrig', 'reward', 'trEnd', 'encodeA', 'encodeB', 'lick')
end

%% task event detection
[trStartIdx,~,~] = detecteventbythreshold(trStart, 25000, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false); % trial Start
[trEndIdx,~,~] = detecteventbythreshold(trEnd, 25000, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1)); % trial End
if length(trStartIdx)==length(trEndIdx)
    if ~unique(trEndIdx - trStartIdx>0)
        error('Trial End and Start indices do not make sense!')
    end
elseif length(trStartIdx)-length(trEndIdx)==1
    if ~unique(trEndIdx - trStartIdx(1,1:length(trEndIdx))>0)
        error('Trial End and Start indices do not make sense!')
    end
else
    error('Trial End and Start indices do not make sense!')
end

rwdIdx     = detecteventbythreshold(reward, 25000, 50, 'stdFactor',1, 'plotRez',false, 'chunkPulses', false, 'detectLater', trStartIdx(1));  % reward
lickIdx    = detecteventbythreshold(lick, 25000, 30, 'stdFactor',3, 'plotRez',false, 'chunkPulses', false);    % lick

[camTrigRiseIdx, camTrigFallIdx, camPulseTrainIdx] = detecteventbythreshold(camTrig, 25000, 2, 'stdFactor', 1, 'plotRez',false, 'chunkPulses', true, 'chunkInterval', 2000); % camera trigger
%camTrigRiseIdx1ms = round(camTrigRiseIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms
%camTrigFallIdx1ms = round(camTrigFallIdx./round(nSamp/1000)); % adjust the time resolution to be 1ms

%% spot the trial-by-trial and all-trials behavioral csv files
if isempty(dir(fullfile(p.Results.filePath,'201*')))
    error('Cannot find the trial-by-trial behavior data csv files!')
end

behFilePath = dir(fullfile(p.Results.filePath,'201*')); % dir where the trial-by-trial behavioral csv files are saved
tbytCsvList = dir(fullfile(behFilePath.folder,behFilePath.name,'trial_*'));    % trial-by-trial files
allTrialCsv = dir(fullfile(behFilePath.folder,behFilePath.name,'trials.csv')); % all trial file
if length(allTrialCsvList)==1
    trialsFileName = fullfile(allTrialCsv.folder,allTrialCsv.name); 
    trialsCsv = readtable(trialsFileName);
else
    error('More than one trials.csv file detected!')
end

[~,tbytCsvdateSort] = sort(datenum({tbytCsvList(:).date}, 'dd-mmm-yyyy hh:MM:ss'), 1, 'ascend'); % sorted fileList

%% stepper encoder data; pin A, pin B
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

for t = 1:length(trStartIdx) % increment trials
    if ~isempty(find(trEndIdx>trStartIdx(t),1)) && trEndIdx(find(trEndIdx>trStartIdx(t),1))>trStartIdx(t) % if there's a trEnd
        % redefine the trial start as the joystick in position timepoint (trJsReadyTime) by examining the baseline encoder data,
        % as the trStart defined with the cue tone onset doesn't align
        % perfectly with the joystick being in the start position.
        jsRez(t).trStart = trStartIdx(t); % trStart detected by the go-cue onset
        trBaseRange = trStartIdx(t)-2*nSamp:trStartIdx(t)+(nSamp/10/2)-1; % 2 sec baseline, just take 2 sec before the cue onset, with padding on the righthand side corresponding to 50 ms
        trBasePinA = encodeA(trBaseRange); % pin A this trial baseline
        trBasePinB = encodeB(trBaseRange); % pin B this trial baseline
        biTrBasePinA = double(trBasePinA>pinAmidPoint); % binarize the pinA signal
        biTrBasePinB = double(trBasePinB>pinBmidPoint); % binarize the pinB signal
        
        [trBaseJsState, ~] = readStepperEncoder(biTrBasePinA, biTrBasePinB); % read out the Js position from the binarized quadrature encoded signals
        [trBaseLatestStillTime] = findJsReadyPt(cumsum(trBaseJsState), 100, nSamp); % find the Js ready time (trStart defined as the cue onset doesn't match with the moment when the Js is in position)
        
        if isnan(trBaseLatestStillTime)    
        elseif trStartIdx(t)-2*nSamp+trBaseLatestStillTime <= jsRez(t).trStart
            jsRez(t).trJsReadyTime = trStartIdx(t)-2*nSamp+trBaseLatestStillTime; % get the trJsReady time
            jsRez(t).trEnd = trEndIdx(find(trEndIdx > trStartIdx(t),1,'first'));  % get the trEnd time
            
            trRange = jsRez(t).trJsReadyTime:jsRez(t).trEnd; % trial range aligned to the Js ready time
            
            trPinA = encodeA(trRange); % pin A this trial aligned to the Js ready
            trPinB = encodeB(trRange); % pin A this trial aligned to the Js ready
            
            % binarize the encoder data
            biTrPinA = double(trPinA>pinAmidPoint); % binarize the encoder signal
            biTrPinB = double(trPinB>pinBmidPoint); % binarize the encoder signal 
            [trJsState, ~] = readStepperEncoder(biTrPinA, biTrPinB); % read out the Js position from the binarized quadrature encoded signals
            jsRez(t).trJsTraj = cumsum(trJsState);  % get the cumulative joystick trajectory for this trial
            jsRez(t).dctrJsTraj = decimate(jsRez(t).trJsTraj ,round(nSamp/1000)); % decimate the Traj (downsampling the 25kHz data at 1kHz, so the time resoluation to be 1ms)
            
            % filter the decimated Js Traj
            sgfiltFramelen = floor(length(jsRez(t).dctrJsTraj)/(p.Results.trialTimeout)*p.Results.sgfiltFramelen);
            if mod(sgfiltFramelen,2)==0 % the frame length for sgolayfilt needs to be an odd positive integer
                sgfiltFramelen = sgfiltFramelen+1;
            end
            
            if sgfiltFramelen <= 3 % the order of polynomial fit for the sgolayfilt needs to be less than the frame length
                sgfiltFramelen = 5;
            end
            jsRez(t).smdctrJsTraj = sgolayfilt(jsRez(t).dctrJsTraj,3,sgfiltFramelen);
            
            % read trial-by-trial behavioral csv files
            tempFileName = fullfile(behFilePath.folder,behFilePath.name,tbytCsvList(tbytCsvdateSort(t)).name);
            tempTimeTraj = csvread(tempFileName,1,1);  % csv read with row offset (header) and column offset (date)
            jsRez(t).csvTime = tempTimeTraj(:,1); % csv time in ms
            jsRez(t).csvTraj = tempTimeTraj(:,2); % csv joystick trajectory
            jsRez(t).csvItpTraj = interp1(1:length(jsRez(t).csvTraj),jsRez(t).csvTraj,linspace(1,length(jsRez(t).csvTraj),length(jsRez(t).smdctrJsTraj)));
            
            % determine if the trial got rewarded
            if ~isempty(find(abs(rwdIdx-jsRez(t).trEnd)<=nSamp,1)) % in case there's reward delivery within 1-sec window relative to the trial end
                jsRez(t).rewarded = true;
            else
                jsRez(t).rewarded = false;
            end
            
            % classfy the trial ('sp': successfull pull, 'ps': push, 'im': immature pull, 'to': timeout, 'nn': not identified)
            if jsRez(t).rewarded % if rewarded
                if ~isempty(find(jsRez(t).smdctrJsTraj<trialsCsv.pull_threshold(t),1)) % check the negative threshold crossing
                    jsRez(t).trialType = 'sp';
                end
            else % if not rewarded
                if length(jsRez(t).smdctrJsTraj)>p.Results.trialTimeout-100 % timeout
                    jsRez(t).trialType = 'to';
                elseif isempty(find(jsRez(t).smdctrJsTraj<trialsCsv.pull_threshold(t),1)) && ~isempty(find(jsRez(t).smdctrJsTraj>20,1)) % push (no pull beyond the pull threshold && push beyond a certain threshold)
                    jsRez(t).trialType = 'ps';
                elseif ~isempty(find(jsRez(t).smdctrJsTraj<trialsCsv.pull_threshold(t),1)) % pull (unrewarded)
                    if ~isempty(find(jsRez(t).smdctrJsTraj>20,1))
                        jsRez(t).trialType = 'impullandpush';
                    else
                        jsRez(t).trialType = 'im';
                    end
                else
                    jsRez(t).trialType = 'nn';
                end
            end
            
            
                
            
            
            
            
            
            %figure; plot(interp1(1:length(jsRez(t).csvTraj),jsRez(t).csvTraj,1:length(jsRez(t).dcsmtrJsTraj)));
        end
    else % if there's no trialEnd left
    end
    fprintf('completed trial #%d\n', t);
end

figure; hold on; plot(jsRez(t).smdctrJsTraj); 
figure; plot(diff([0 jsRez(t).smdctrJsTraj]));
figure; plot(diff([0 diff([0 jsRez(t).smdctrJsTraj])]));
plot(jsRez(125).csvItpTraj)

%% Position/velocity data


% Save relevant BehVariables
cd(filePath)
if p.Results.laserUsed
    save('BehVariables','Xpos','Ypos','positionData','lick','sole','reach0','pos1','pos2','xpos1','ypos1','xpos2','ypos2','vel1','vel2','ts','p', 'laser', 'pseudoLaser' ) % append the position/velocity data variables
else
    save('BehVariables','Xpos','Ypos','positionData','lick','sole','reach0','pos1','pos2','xpos1','ypos1','xpos2','ypos2','vel1','vel2','ts','p' ) % append the position/velocity data variables
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_Js( filePath, vargs )
        % parse input, and extract name-value pairs
        default_numbNeuralProbe = 0;  % specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
        default_numbChEachProbe = 64; % specify how many channels are on the probe
        default_trStartCh = 33; % ch# for trial start
        default_camTrigCh = 34; % ch# for camera trigger
        default_rewardCh  = 35; % ch# for reward delivery
        default_trEndCh   = 36; % ch# for trial end
        default_encodeACh = 37; % ch# for stepper encoder A
        default_encodeBCh = 39; % ch# for stepper encoder B
        default_lickCh    = 1;  % ch# for lick detect (unattenuated channel)
        default_sgfiltFramelen = 101; % frame length for the sgolayfilt
        default_trialTimeout = 10000; % trial timeout duration
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'numbNeuralProbe',default_numbNeuralProbe)
        addParameter(p,'numbChEachProbe',default_numbChEachProbe)
        addParameter(p,'trStartCh',default_trStartCh)
        addParameter(p,'camTrigCh',default_camTrigCh)
        addParameter(p,'rewardCh',default_rewardCh)
        addParameter(p,'trEndCh',default_trEndCh)
        addParameter(p,'encodeACh',default_encodeACh)
        addParameter(p,'encodeBCh',default_encodeBCh)
        addParameter(p,'lickCh',default_lickCh)
        addParameter(p,'sgfiltFramelen',default_sgfiltFramelen)
        addParameter(p,'trialTimeout',default_trialTimeout)
        
        parse(p,filePath,vargs{:})
    end


end

