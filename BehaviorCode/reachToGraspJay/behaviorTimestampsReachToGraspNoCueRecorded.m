function behaviorTimestampsReachToGraspNoCueRecorded(filePath)
cd(filePath)

binFile = dir(fullfile(filePath,'*.nidq.bin')); % look for nidq.bin file

if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Parse the corresponding metafile
meta = ReadMeta(binName, filePath); % get the meta data (structure)

% Read the binary data (entire samples)
nSamp         = SampRate(meta);          % sampling rate (default: 25kHz)
totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds

% Specify the relevant behavioral channel numbers
camTrigCh = 64+1;
cueCh = 64+2;
wvsCh = 64+3;
laserCh = 64+4;

if exist('gaincorrectRawData.mat','file')==2
    load(fullfile(filePath,'gaincorrectRawData.mat'))
else
    % preallocate the behavioral data arrays
    camTrig = zeros(1,floor(totalTimeSecs*1000)); % the time resolution will be 1000Hz (1ms) after decimation
    cue = zeros(1,floor(totalTimeSecs*1000)); %
    wvs = zeros(1,floor(totalTimeSecs*1000)); %
    laser = zeros(1,floor(totalTimeSecs*1000)); %
    
    for i = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
        tempDataArray = ReadBin(i*nSamp, nSamp, meta, binName, filePath); % read bin data for each second
        tempCamTrig  = decimate(tempDataArray(camTrigCh,:),round(nSamp/1000)); % decimate the data
        tempCue  = decimate(tempDataArray(cueCh,:),round(nSamp/1000));
        tempWvs  = decimate(tempDataArray(wvsCh,:),round(nSamp/1000));
        tempLaser = decimate(tempDataArray(laserCh,:),round(nSamp/1000));
        
        camTrig(1,i*1000+1:(i+1)*1000) = tempCamTrig; % camera trigger (concatenated the decimated data)
        cue(1,i*1000+1:(i+1)*1000) = tempCue; %  for trial onset
        wvs(1,i*1000+1:(i+1)*1000) = tempWvs; % wave surfer signal
        laser(1,i*1000+1:(i+1)*1000) = tempLaser; % laser (1. cue+laser(2s), 2. laser-only (2s), 3. tagg-laser (1s))
        
        clearvars temp*
        fprintf('processed %d\n', i+1)
    end
    clearvars i
    
    % Gain correction for channnels of interest
    camTrig = GainCorrectNI(camTrig, 1, meta);   % gain-corrected voltage trace for camTrig
    cue = GainCorrectNI(cue, 1, meta);   % gain-corrected voltage trace for cue
    wvs = GainCorrectNI(wvs, 1, meta);   % gain-corrected voltage trace for wvs
    laser = GainCorrectNI(laser, 1, meta);   % gain-corrected voltage trace for laser
    
    save(fullfile(filePath,'gaincorrectRawData'),'camTrig','cue','wvs','laser')
end

%% camera trigger detection
[camRiseIdx, camFallIdx, ~] = detecteventbythreshold(camTrig, 1000, 4000, 'stdFactor', 1, 'plotRez',false, 'chunkPulses', false, 'chunkInterval', 1000, 'correctLongPulse', false); % camera trigger

%% laser detection
% detect rise and fall of laser
[laserRiseIdx, laserFallIdx, ~] = detecteventbythreshold(laser, 1000, 5, 'stdFactor', 2, 'plotRez',false, 'chunkPulses', false, 'chunkInterval', 1000, 'correctLongPulse', false); % camera trigger

if length(laserRiseIdx)~=length(laserFallIdx)
    error('The lengths of laser on and off differ!')
end

laserDur = laserFallIdx-laserRiseIdx;
%if sum(laserDur<500)>0
%    error('Check out the laser duration!')
%end
%laserDur = laserDur(laserDur>500);

if laserRiseIdx(1)-camRiseIdx(3)<2000
    laserOnlyIdx = false(1,length(laserRiseIdx)); 
    laserOnlyIdx(1:2:length(laserRiseIdx)) = true; 
    laserCueIdx = ~laserOnlyIdx; 
end

ts.tagLaser1s = laserRiseIdx(laserDur>500&laserDur<1100&laserOnlyIdx);
%laser2s = laserRiseIdx(laserDur>1900);
%laserRiseInt = arrayfun(@(x) x-1000:x+1000, laserRiseIdx(laserDur>1900), 'un', 0);
%laserCueI = cell2mat(cellfun(@(a) sum(ismember(valCueIdx,a))==1, laserRiseInt, 'un', 0));
ts.laserCue2s = laserRiseIdx(laserDur>1900&laserCueIdx); %laser2s(laserCueI);
ts.laserOnly2s = laserRiseIdx(laserDur>1900&laserOnlyIdx); %laser2s(~laserCueI);
fprintf('completed laser detection!');

%% cue detection (WHEN NO CUE WAS RECORDED BY MISTAKE, TRY TO INFER FROM CAMTRIG AND LASER)
% [~,cueStd,~] = meanstdsem(abs(cue)');
% cueThres      = mean(abs(cue))+2*cueStd;
% %cueThres     = 0.8*10^-4; 
% cueIdx        = find(abs(cue)>cueThres);
% valCueIdx     = cueIdx(diff([0,cueIdx])>5000);
% ts.cue = valCueIdx;
% fprintf('Cues detected: %d\n', length(valCueIdx));

cueIdx = camRiseIdx+1000; 
valCueIdx = cueIdx(diff([0,cueIdx])>5000);
ts.cue = valCueIdx;
fprintf('Cues detected: %d\n', length(valCueIdx));

% if length(camRiseIdx)==length(camFallIdx)
%     if sum(camFallIdx-camRiseIdx>0)==length(camRiseIdx)
%         camTrigC = mat2cell([camRiseIdx',camFallIdx'],ones(length(camRiseIdx),1),2); 
%         camTrigC = cellfun(@(a) a(1):a(end),camTrigC,'un',0); 
%         laserOnlyCamI = cell2mat(cellfun(@(a) sum(ismember(laserRiseIdx(laserOnlyIdx),a))>0, camTrigC,'un',0));
%     end
% end
%ts.cue = camRiseIdx(~laserOnlyCamI)+1000; 

% Save relevant BehVariables
cd(filePath)
save(fullfile(filePath,'BehVariablesJayRtg'),'ts','camTrig','laser','cue','wvs') % append the position/velocity data variables

end

