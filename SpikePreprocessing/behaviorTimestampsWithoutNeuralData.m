function behaviorTimestampsWithoutNeuralData(filePath,varargin)
%behaviorTimestamps takes filePath in which the nidq.bin file exists, and using
% the nidq.bin and nidq.meta files, occurrences of a set of behavioral events - 
% reward delivery, joystick movement trajectories, licks, laser delivery
% are detected and stored in the filePath.
% Inputs to the function can be provided as name-value pair arguments.
% e.g. Run, 
% behaviorTimestamps(filePath,'numbNeuralProbe',0,'numbChEachProbe',64,...
% ...,'XposCh',33,'YposCh',37,'soleCh',3,'lickCh',5,'laserCh',7,'numbTagLasers',30,'artifactRmv',true)
% modified on 5/17/18 to process the pseudoLaser signal (pseudoLaser pulses without actual laser delivery)

%addpath(genpath('/Volumes/RAID2/parkj/MATLAB'))
p = parse_input_beh(filePath, varargin ); % parse input
%p = parse_input_beh(filePath, {} ); % use this line instead when running line-by-line 

cd(p.Results.filePath)

if ~isempty(dir('BehVariables.mat')) % if the BehVariables.mat file already exists in the filePath
    answer = questdlg('BehVariable.mat already exists, Would you like to replace it?','Choice','Replace','Cancel','Cancel');
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
meta = ReadMeta(binName, p.Results.filePath); % get the meta data (structure)

% Read the binary data (entire samples)
nSamp         = SampRate(meta);          % sampling rate (default: 25kHz) 
totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds

% Specify the relevant behavioral channel numbers
%NumbNeuralProbe = 0; % Specify how many probes were used (e.g. zero if no NI neural probe was used)
XposCh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.XposCh;  % channel # for X position (default channel numbers for 64 channel recording)
YposCh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.YposCh;  % channel # for Y position
soleCh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.soleCh;  % channel # for solenoid (water reward delivery) 
lickCh  = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.lickCh;  % channel # for lick port
laserCh = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.laserCh; % channel # for laser (laser TTL) 
pseudoLaserCh = p.Results.numbChEachProbe*p.Results.numbNeuralProbe+p.Results.pseudoLaserCh; % channel # for pseudoLaser (laser TTL triggerred by joystick displacement crossing the laser stim threshold with and without the actual laser delivery) 

% preallocate the behavioral data arrays
Xpos    = zeros(1,floor(totalTimeSecs*1000)); % the time resolution will be 1000Hz (1ms) after decimation
Ypos    = zeros(1,floor(totalTimeSecs*1000)); %
lick    = zeros(1,floor(totalTimeSecs*1000)); %
sole    = zeros(1,floor(totalTimeSecs*1000)); %
laser   = zeros(1,floor(totalTimeSecs*1000)); %
pseudoLaser = zeros(1,floor(totalTimeSecs*1000)); %

for i = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
    tempDataArray = ReadBin(i*nSamp, nSamp, meta, binName, p.Results.filePath); % read bin data for each second
    tempXpos  = decimate(tempDataArray(XposCh,:),round(nSamp/1000)); % decimate the data
    tempYpos  = decimate(tempDataArray(YposCh,:),round(nSamp/1000));
    templick  = decimate(tempDataArray(lickCh,:),round(nSamp/1000));
    tempsole  = decimate(tempDataArray(soleCh,:),round(nSamp/1000));
    templaser = decimate(tempDataArray(laserCh,:),round(nSamp/1000));
    temppseudoLaser = decimate(tempDataArray(pseudoLaserCh,:),round(nSamp/1000));
    Xpos(1,i*1000+1:(i+1)*1000) = tempXpos; % accumulated the decimated data
    Ypos(1,i*1000+1:(i+1)*1000) = tempYpos;
    lick(1,i*1000+1:(i+1)*1000) = templick;
    sole(1,i*1000+1:(i+1)*1000) = tempsole;
    laser(1,i*1000+1:(i+1)*1000) = templaser;
    pseudoLaser(1,i*1000+1:(i+1)*1000) = temppseudoLaser;
    clearvars temp*
    fprintf('processed %d\n', i+1) 
end
clearvars i

% Gain correction for channnels of interest
if strcmp(meta.typeThis, 'imec') % in case recording via imec 
    Xpos = GainCorrectIM(Xpos, 1, meta);   % gain-corrected voltage trace for Xpos 
    Ypos = GainCorrectIM(Ypos, 1, meta);   % gain-corrected voltage trace for Ypos 
    lick = GainCorrectIM(lick, 1, meta);   % gain-corrected voltage trace for lick
    sole = GainCorrectIM(sole, 1, meta);   % gain-corrected voltage trace for solenoid
    laser = GainCorrectIM(laser, 1, meta); % gain-corrected voltage trace for laser
    pseudoLaser = GainCorrectIM(pseudoLaser, 1, meta); % gain-corrected voltage trace for pseudolaser
else    % in case of recording via NI board
    Xpos = GainCorrectNI(Xpos, 1, meta);   % gain-corrected voltage trace for Xpos 
    Ypos = GainCorrectNI(Ypos, 1, meta);   % gain-corrected voltage trace for Ypos 
    lick = GainCorrectNI(lick, 1, meta);   % gain-corrected voltage trace for lick 
    sole = GainCorrectNI(sole, 1, meta);   % gain-corrected voltage trace for solenoid
    laser = GainCorrectNI(laser, 1, meta); % gain-corrected voltage trace for laser
    pseudoLaser = GainCorrectNI(pseudoLaser, 1, meta); % gain-corrected voltage trace for pseudolaser
end    

%% reward (digital pulses for solenoid activation) detection 
[~,soleStd,~] = meanstdsem(abs(sole)');       % std of the lick input signal
rewThres      = mean(abs(sole))+soleStd;      % this seems to work as a reasonable threshold for detecting lick stim
rewIdx        = find(abs(sole)>rewThres);     % find points crossing the lick threshold
valRewIdx     = rewIdx(diff([0,rewIdx])>500);  % this prevents redundant detections
% hold on; plot(sole); plot(valRewIdx,rewThres,'or'); hold off

%% Position/velocity data
if p.Results.artifactRmv % in case artifact remove is true
    denoiseXpos=denoiseByTemplateSubtraction(Xpos,valRewIdx);
    denoiseYpos=denoiseByTemplateSubtraction(Ypos,valRewIdx);
    pseudoLaser=denoiseByTemplateSubtraction(pseudoLaser,valRewIdx); 
    positionData = [denoiseXpos; denoiseYpos]; % denoised (without the solenoid artifact) joystick position data
else % in case artifact remove is not selected
    positionData = [Xpos; Ypos]; % joystick position data
end

if p.Results.reachBeforeLastReward
    positionData = positionData(:,1:valRewIdx(end));
else
end

% get reach properties
[ reachStart, reachStop, reach0, pos1, pos2, xpos1, ypos1, xpos2, ypos2 ] = getReachTimesJP( positionData );     % all reach traces, aligned to start (pos1), to stop (pos2)

vel1 = diff(pos1, 1, 2); % reach velocity aligned to reach start(differentiation of pos1) 
vel2 = diff(pos2, 1, 2); % reach velocity aligned to reach stop (differentiation of pos2)
% plot(reach0) % reachMW is the amplitude readout of the whole session

%% Get other task events - reward delivery, licks, laser stimulation
% lick (digital pulses for lick) detection 
[~,lickStd,~] = meanstdsem(abs(lick)');          % std of the lick input signal
lickThres     = mean(abs(lick))+2.5*lickStd;         % this seems to work as a reasonable threshold for detecting lick stim
lickIdx       = find(abs(lick)>lickThres);       % find points crossing the lick threshold
valLickIdx    = lickIdx(diff([0,lickIdx])>10);   % this prevents redundant detections

%this part visualizes the detected lick traces for validation
valLickCount = 0;
for i = 1:length(valLickIdx)
    %hold on;
    if valLickIdx(i)-100>0 && valLickIdx(i)+100<length(lick)
        valLickCount = valLickCount + 1;
        lickTraces(valLickCount,:) = lick(1,valLickIdx(i)-100:valLickIdx(i)+100);
        %plot(valLickIdx(i)-100:valLickIdx(i)+100,lickTraces(valLickCount,:),'m')
        %plot(valLickIdx(i), lick(valLickIdx(i)),'c*')
    else
        
    end
end
clearvars i
% validation with plot
%hold on; plot(lick); plot(valLickIdx,lick(valLickIdx),'or'); hold off

if p.Results.laserUsed
    % laser (TTL pulses for laser) detection
    [~,laserStd,~] = meanstdsem(abs(laser)');             % std of the laser input signal
    laserThres     = mean(abs(laser))+laserStd;           % this seems to work as a reasonable threshold for detecting laser stim
    laserIdx       = find(abs(laser)>laserThres);         % find points crossing the laser threshold
    valLaserIdx    = laserIdx(diff([0,laserIdx])>1000);   % this prevents redundant detections
    stmLaser       = valLaserIdx;  % randomly selected reach-evoked opto-stimulations
    
    if p.Results.reachBeforeLastReward
        stmLaser = stmLaser(stmLaser<valRewIdx(end)); % only take stmLaser occurred before the last reward delivery
    else
    end
    
    % validation with plot
    %hold on; plot(laser); plot(stmLaser,laserThres,'or'); plot(tagLaser,laserThres,'og'); hold off;
    
    % pseudoLaser detection
    rezeroPseudoLaser = pseudoLaser - mean(pseudoLaser); % re-zero pseudoLaser It's a nice practice to rezero the digital signals
    [~,rezeroPseudoLaserStd,~] = meanstdsem(abs(rezeroPseudoLaser)');         % std of the laser input signal
    pseudoLaserThres     = mean(abs(rezeroPseudoLaser))+rezeroPseudoLaserStd; % this seems to work as a reasonable threshold for detecting laser stim
    pseudoLaserIdx       = find(abs(rezeroPseudoLaser)>pseudoLaserThres);     % find points crossing the laser threshold
    valpseudoLaserIdx    = pseudoLaserIdx(diff([0,pseudoLaserIdx])>2000);     % this prevents redundant detections
    
    % detect reaches with stim on (also the simulations delivered during reach or not)
    stimReachIdx = zeros(length(reachStart),1);     % index for stim on or not to be used for reachStart or reachStop
    stimReachLaserIdx = zeros(length(stmLaser),1);  % index for stmLaser to be delivered during completed reach or not
    
    for i = 1:length(reachStart)
        
        tmpStmLaserId = find(stmLaser > reachStart(i),1,'first'); % the first stim on after the current reach start
        if stmLaser(tmpStmLaserId) <= reachStop(i)
            stimReachIdx(i) = true; % mark it as a stimulated trial
            stimReachLaserIdx(tmpStmLaserId) = true; % mark it as a stimulation delivered during a complete reach
        end
    end
end

% Build a structure for timestamps
ts.reachStart = reachStart;     % reachStart
ts.reachStop  = reachStop;      % reachStop
ts.reward     = valRewIdx;      % reward deliveries
ts.lick       = valLickIdx;     % licks

if p.Results.laserUsed
    ts.laser      = valLaserIdx;    % laser stimulations all lasers
    ts.stmLaser   = stmLaser;       % randomly selected reach-evoked opto-stimulations
    ts.stmLaserReach = stmLaser(logical(stimReachLaserIdx)); % stim trials occurred during full reaches
    ts.pseudoLaser = valpseudoLaserIdx; % pseudoLaser TTL pulses (without actual laser deliveries)
    ts.stmReachStart = reachStart(logical(stimReachIdx)); % reachStart of stim on trials
    ts.stmReachStop  = reachStop(logical(stimReachIdx));  % reachStop of stim on trials
end

% Reach rewarded or not
ts.reachRew = zeros(length(ts.reachStart),1); % 1st column: reward logic
for t = 1:length(ts.reachStart)
    if ~isempty(find(ts.reward > ts.reachStart(t),1)) % in case, there's any reward delivery after this reach
        nextRew = ts.reward(min(find(ts.reward > ts.reachStart(t),1))); % the very next reward after this reachStart
        if t<length(ts.reachStart)
            ts.reachRew(t,1) = nextRew < ts.reachStart(t+1); % see if a reward was delivered before the next reachStart or not
        elseif t==length(ts.reachStart)
            ts.reachRew(t,1) = true;
        end
    else % in case there's no subsequent reward delivery at all  
        ts.reachRew(t,1) = false;
    end
end

% Save relevant BehVariables
cd(filePath)

if p.Results.laserUsed
    save('BehVariables','Xpos','Ypos','positionData','lick','sole','lickTraces','reach0','pos1','pos2','xpos1','ypos1','xpos2','ypos2','vel1','vel2','ts','p', 'laser', 'pseudoLaser' ) % append the position/velocity data variables
else
    save('BehVariables','Xpos','Ypos','positionData','lick','sole','lickTraces','reach0','pos1','pos2','xpos1','ypos1','xpos2','ypos2','vel1','vel2','ts','p' ) % append the position/velocity data variables
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_beh( filePath, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function 
        % parse input, and extract name-value pairs
                
        default_numbNeuralProbe = 0;  % Specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
        default_numbChEachProbe = 64; % Specify the number of sites on the NIboard probe
        default_XposCh = 33; % channel # for X position (default channel numbers for 64 channel recording)
        default_YposCh = 37; % channel # for Y position
        default_pseudoLaserCh = 39; % channel # for pseudoLaser Pulses (without actual laser delivery)
        default_soleCh = 3;  % channel # for solenoid (water reward delivery)
        default_lickCh = 5;  % channel # for lick port
        default_laserCh = 7; % channel # for laser (laser TTL)
        default_numbTagLasers = 30; % the number of tagging trials given at the end of the experiment
        default_artifactRmv = true; % if true, removes the solenoid artifact from Xpos and Ypos channels by template subtraction
        default_reachBeforeLastReward = true; % logical to detect reaches before the last reward delivery
        default_laserUsed = true; % logical to indicate whether laser stimulation was used during the session or not
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'numbNeuralProbe',default_numbNeuralProbe)
        addParameter(p,'numbChEachProbe',default_numbChEachProbe)
        addParameter(p,'XposCh',default_XposCh)
        addParameter(p,'YposCh',default_YposCh)
        addParameter(p,'pseudoLaserCh', default_pseudoLaserCh)
        addParameter(p,'soleCh',default_soleCh)
        addParameter(p,'lickCh',default_lickCh)
        addParameter(p,'laserCh',default_laserCh)
        addParameter(p,'numbTagLasers',default_numbTagLasers)
        addParameter(p,'artifactRmv',default_artifactRmv)
        addParameter(p,'reachBeforeLastReward',default_reachBeforeLastReward)
        addParameter(p,'laserUsed', default_laserUsed)
        
        parse(p,filePath,vargs{:})
        
    end

end

