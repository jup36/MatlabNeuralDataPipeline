function js2p0_tbytSpkHandJsPreprocess_50ms_stimPstim(filePath)
%This is a preprocessing function to first demarcate trials into blocks of
% different joystick load & position combinations using the function 'jkvtBlockParse'.
% Then it gets trial-by-trial binned (e.g. 20-ms bin) spike count matrices aligned to
% an event ('reachStart' or 'pullStart' or 'trJoystickReady') depending on the trial type.
% For each trial, hand/joystick trajectories (from DLC) and joystick kinetics
% (force applied from joystick encoder) are acquired, if available, with the same alignment
% as the spike data but curtailed as those kinematic informations get cut
% off right after the end of execution. Additionally, some useful
% trial-by-trial behavioral variables are computed for future decoding
% anlaysis such as reach angle, max force, initial hand position etc.
% The output structure 'ss', 'jkvt', 'trI' are saved in the same filePath.
% filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';

%% 1. Load data
%clc; clearvars; close all;
%filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles';
cd(filePath)
% neural and behavioral data
spkDir_CtxStr = dir('binSpkCountSTRCTX*');
spkDir_Cg = dir('binSpkCountCg*');

if ~isempty(spkDir_CtxStr)
    load(fullfile(spkDir_CtxStr(1).folder, spkDir_CtxStr(1).name),'spkTimesCell','jkvt')
    spkTimesCell_CtxStr = spkTimesCell; clearvars spkTimesCell
    % cortex/striatum index
    depth = cell2mat(cellfun(@(a) a(2), spkTimesCell_CtxStr(4,:),'un',0))'; % depth from pial surface
    ctxI = depth<1900; % cortex index
    spkTimesCellCTX = spkTimesCell_CtxStr(:,ctxI);
    depthCtx = depth(ctxI);
    strI = depth>2100; % striatum index
    spkTimesCellSTR = spkTimesCell_CtxStr(:,strI);
    depthStr = depth(strI);
end

if ~isempty(spkDir_Cg)
    load(fullfile(spkDir_Cg(1).folder, spkDir_Cg(1).name),'spkTimesCell')
    if ~exist('jkvt', 'var')
        load(fullfile(spkDir_Cg(1).folder, spkDir_Cg(1).name), 'jkvt')
    end

    spkTimesCellCg = spkTimesCell; clearvars spkTimesCell
    depthCg = cell2mat(cellfun(@(a) a(2), spkTimesCellCg(4,:),'un',0))'; % depth from pial surface
end

%% align hand trajectories to neural data
vfT = {jkvt(:).vFrameTime}'; % video frame time
hTrj = {jkvt(:).hTrjF}'; % hand trajectory

sm_kernel = TNC_CreateGaussian(250,30,500,20); % a kernel for smoothing (mu, sigma, time, dT)
binSize = 50; % 50 ms
spkBin = -1000:2000;

%% parse trial/block information
% pull torque & reach position combinations
pTqs = round(unique([jkvt(:).pull_torque])./10); % pull torque list (divided by 10 for ordering)
rPos = round(unique([jkvt(:).reachP1])*10); % reach position list (multiplied by 10 for ordering)
posTq = unique(repmat(pTqs,[length(rPos) 1])+repmat(rPos',[1 length(pTqs)])); % position-torque combinations
posTqCnt = zeros(length(posTq),1); % count occurrence of each combination

%% parse block transition and block type
jkvt = jkvtBlockParse(jkvt);

%% align behavioral and neural trajectories
% determine where to align based on trial types
% Trial type index
trI.spI = cell2mat(cellfun(@(a) strcmpi(a,'sp'),{jkvt(:).trialType},'un',0)); % successful pull
trI.toI = cell2mat(cellfun(@(a) strcmpi(a,'to'),{jkvt(:).trialType},'un',0)); % time out
trI.psI = cell2mat(cellfun(@(a) strcmpi(a,'ps'),{jkvt(:).trialType},'un',0)); % push
trI.pmppI = cell2mat(cellfun(@(a) strcmpi(a,'pmpp'),{jkvt(:).trialType},'un',0)); % premature pull and push
trI.pmI = cell2mat(cellfun(@(a) strcmpi(a,'pm'),{jkvt(:).trialType},'un',0)); % premature pull and push

for t = 1:size(jkvt,2)
    %% get a time point to align depending on the trial type
    ss(t).trialType = jkvt(t).trialType;
    % get time point to which neural data and hand trajectories are aligned
    if strcmpi(jkvt(t).trialType, 'sp') % a successful pull trial
        if ~isempty(jkvt(t).rStartToPull) % if there's rStartToPull
            ss(t).timeAlign = jkvt(t).rStartToPull; % align to reach start to pull
            ss(t).evtAlign  = 'rStart';
            if ~isempty(jkvt(t).rStopToPull)
                ss(t).rEnd = jkvt(t).rStopToPull; % reachStop
            elseif ~isempty(jkvt(t).pullStops)
                ss(t).rEnd = jkvt(t).pullStops;   % instead use pullStop
            else
                ss(t).rEnd = jkvt(t).trEnd;
            end
        elseif ~isempty(jkvt(t).pullStarts) && ~isempty(jkvt(t).pullStops)
            ss(t).timeAlign = jkvt(t).pullStarts;
            ss(t).rEnd = jkvt(t).trJsReady+jkvt(t).movKins.pullStop; % reachStop
            ss(t).evtAlign  = 'pStart';
        else
            ss(t).timeAlign = jkvt(t).trJsReady; % if no reachStart detected, just align to the joystick ready
            ss(t).evtAlign  = 'trJsReady';
            ss(t).rEnd = jkvt(t).trEnd;
        end
    else % not a success trial
        if ~isempty(jkvt(t).hTrjRstart) && ~isempty(jkvt(t).hTrjRstop) % if there's a detected reachStart align to that
            ss(t).timeAlign = jkvt(t).vFrameTime(jkvt(t).hTrjRstart(end));
            ss(t).evtAlign  = 'rStart';
            ss(t).rEnd = jkvt(t).vFrameTime(min(jkvt(t).hTrjRstop(end),length(jkvt(t).vFrameTime))); % reachStop
        else
            ss(t).timeAlign = jkvt(t).trJsReady; % if no reachStart detected, just align to the joystick ready
            ss(t).evtAlign  = 'trJsReady';
            ss(t).rEnd = jkvt(t).trEnd;
        end
    end

    % get spike time bins and binned spike count matrices
    ss(t).spkTimeBins = ss(t).timeAlign + (spkBin(1):binSize:spkBin(end-1));
    %spkTime1msBins = ss(t).timeAlign + (spkBin(1):spkBin(end-1));

    if exist('spkTimesCellCTX')==1
        ss(t).unitTimeBCtx = psthBINcellPerTrial( spkTimesCellCTX, ss(t).timeAlign, binSize, [abs(spkBin(1)) abs(spkBin(end))]); % binned spikeCounts aligned to this trial
    end
    if exist('spkTimesCellSTR')==1
        ss(t).unitTimeBStr = psthBINcellPerTrial( spkTimesCellSTR, ss(t).timeAlign, binSize, [abs(spkBin(1)) abs(spkBin(end))]); % binned spikeCounts aligned to this trial
    end
    if exist('spkTimesCellCg')==1
        ss(t).unitTimeBCg = psthBINcellPerTrial( spkTimesCellCg, ss(t).timeAlign, binSize, [abs(spkBin(1)) abs(spkBin(end))]); % binned spikeCounts aligned to this trial
    end

    % get unitTimeBin aligned to stimOn (jkvt(t).stimLaserOn) or pseudoStimOn (jkvt(t).pLaserOn)
    if ~isfield(jkvt, 'pLaserOn')
        stimI = ~isnan([jkvt.stimLaserOn]);
        stimOnRelToTrStart = nanmean(cell2mat(cellfun(@(a, b) a-b, {jkvt.trStart}, {jkvt.stimLaserOn}, 'UniformOutput', false))); 
        for jj = find(~stimI)
            jkvt(jj).pLaserOn = jkvt(jj).trStart-stimOnRelToTrStart; 
        end
    end
    
    if ~isnan(jkvt(t).stimLaserOn) % stim trial
        % select trials that had laser on for at least 2 sec
        if jkvt(t).stimLaserOff - jkvt(t).stimLaserOn > 2000
            if exist('spkTimesCellCTX')==1
                ss(t).utbCtxStimAlign = psthBINcellPerTrial( spkTimesCellCTX, jkvt(t).stimLaserOn, binSize, [0 2000]); % binned spikeCounts aligned to this trial
            end
            if exist('spkTimesCellSTR')==1
                ss(t).utbStrStimAlign = psthBINcellPerTrial( spkTimesCellSTR, jkvt(t).stimLaserOn, binSize, [0 2000]); % binned spikeCounts aligned to this trial
            end
            if exist('spkTimesCellCg')==1
                ss(t).utbCgStimAlign = psthBINcellPerTrial( spkTimesCellCg, jkvt(t).stimLaserOn, binSize, [0 2000]); % binned spikeCounts aligned to this trial
            end
        end
    else % control trial
        if isempty(jkvt(t).rStartToPull)
            takePstimTrI = true;
        elseif jkvt(t).rStartToPull-jkvt(t).pLaserOn>2000 % ensure to exclude trials where reach initiated
            takePstimTrI = true;
        end

        if takePstimTrI
            if exist('spkTimesCellCTX')==1
                ss(t).utbCtxPstimAlign = psthBINcellPerTrial( spkTimesCellCTX, jkvt(t).pLaserOn, binSize, [0 2000]); % binned spikeCounts aligned to this trial
            end
            if exist('spkTimesCellSTR')==1
                ss(t).utbStrPstimAlign = psthBINcellPerTrial( spkTimesCellSTR, jkvt(t).pLaserOn, binSize, [0 2000]); % binned spikeCounts aligned to this trial
            end
            if exist('spkTimesCellCg')==1
                ss(t).utbCgPstimAlign = psthBINcellPerTrial( spkTimesCellCg, jkvt(t).pLaserOn, binSize, [0 2000]); % binned spikeCounts aligned to this trial
            end
        end
    end

    %% get interpolated/binned hand position, velocity and force measured from the joystick encoder (all traj aligned to t1n e.g., -1000ms from rStart)
    if ~isempty(hTrj{t}) && ~isempty(ss(t).timeAlign) % if hTrj available
        spikeT = spkBin+ss(t).timeAlign; % 1-ms spike time bins
        [inthTrj, intX] = interpsm(vfT{t},hTrj{t}); % interpolate and smooth
        ss(t).hTrjBfull = inthTrj(:,1:binSize:size(inthTrj,2));
        t1n = ss(t).timeAlign-abs(spkBin(1)); % t1 for neural spike trains
        tEn = ss(t).timeAlign+abs(spkBin(2)); % tE for neural spike trains
        t1h = intX(1); % t1 for hand trajectory
        tEr = ss(t).rEnd;  % time reach ends

        if t1n<t1h % need to extrapolate to the left (earlier)
            extX1 = t1n:intX(end);
            [exthTrj1] = extm(intX,extX1,inthTrj,'linear'); % linear extrapolation
            hTrj1 = exthTrj1(:,1:find(extX1==tEr));
        else
            hTrj1 = inthTrj(:,find(intX==t1n):find(intX==tEr));
        end

        % bin hTrj1
        ss(t).hTrjB = hTrj1(:,1:binSize:size(hTrj1,2)); % bin hand trajectory
        hVel1 = (hTrj1(:,2:end)-hTrj1(:,1:end-1))*(1000/1)/10; % velocity (cm/s)
        ss(t).hVelB = hVel1(:,1:binSize:size(hVel1,2)); % bin hand velocity
        ss(t).hInitPos = nanmedian(ss(t).hTrjB(:,1:max(1,abs(spkBin(1))/binSize/2)),2); % initial hand position
        hDistFromInitPos = cell2mat(cellfun(@(a) sum(sqrt((a-repmat(ss(t).hInitPos,1,size(a,2))).^2),1), {hTrj1}, 'un', 0));  % compute the distance from the jsXY
        ss(t).hDistFromInitPos = hDistFromInitPos(:,1:binSize:size(hTrj1,2)); % bin hand trajectory

        % force trace applied onto the joystick
        if strcmpi(jkvt(t).trialType, 'sp') && isfield(jkvt(t).movKins,'forceMN')
            tmpForce = jkvt(t).movKins.forceMN(1,:); % current trial's smoothed Js force trace (from the encoder, acceleration + Mass)
            tmpForceN = zeros(1, length(tmpForce)); % current trial's force trace only to include negative (means pull) portions
            tJeNI = jkvt(t).movKins.forceMN(1,:)<0; % time index where force applied toward the pull direction
            tmpForceN(tJeNI) = tmpForce(tJeNI);     % pull force trace
            sTmpForceN = conv(tmpForceN,sm_kernel,'same'); % smoothing with convolution
            tJe = jkvt(t).trJsReady:jkvt(t).trEnd-1; % time points for the joystick trace from encoder
            tmpForceNt1ntEr = targetintoreftime(sTmpForceN,tJe,t1n:tEr); % align pulling force to the time frame of hTrj
            ss(t).maxPullForce = min(tmpForceNt1ntEr(1,abs(spkBin(1)):end)); % max pull force (pulls are negative-valued)
            ss(t).forceB = tmpForceNt1ntEr(:,1:binSize:size(hTrj1,2)); % bin the forceTrace
        end
        % get the pullStart and pullStop point and pull index, if exists
        if isfield(jkvt(t).movKins,'pullStart') && isfield(jkvt(t).movKins,'pullStop')
            pullStart = jkvt(t).trJsReady + jkvt(t).movKins.pullStart;
            pullStop = jkvt(t).trJsReady + jkvt(t).movKins.pullStop;
            ss(t).tPullStart = pullStart;
            ss(t).tPullStop = pullStop;
            ss(t).spkPullIdx = pullStart<=ss(t).spkTimeBins & ss(t).spkTimeBins<=pullStop;
            ss(t).spkRchIdx  = ss(t).timeAlign<=ss(t).spkTimeBins & ss(t).spkTimeBins<=ss(t).rEnd;
            ss(t).rchSpeed1ms = sqrt(sum(hVel1(:,max(1,abs(spkBin(1))-100):length(t1n:pullStart)).^2,1)); % cm/s speed (not velocity) during reach phase
            ss(t).maxRchSpeed = max(ss(t).rchSpeed1ms,[],2); % cm/s max speed during reach phase
            ss(t).maxRchSpeedXYZ = max(abs(hVel1(:,max(1,abs(spkBin(1))-100):length(t1n:pullStart))),[],2); % max speed during reach phase on X,Y,Z separately
            hTrjBfullP1_1ms = inthTrj(:,intX<=pullStart);
            ss(t).hTrjBfullP1 = hTrjBfullP1_1ms(:,1:binSize:size(hTrjBfullP1_1ms,2));
        end

        %% get joystick trajectory (all traj aligned to t1n e.g., -1000ms from rStart)
        if ~isempty(jkvt(t).jsTrjValB)
            tJv = vfT{t}(length(vfT{t})-length(jkvt(t).jsTrjValB)+1:end); % joystick traj time points (4ms interval)
            [intjTrj, intjsX] = interpsm(tJv,jkvt(t).jsTrjValB); % interpolate and smooth
            if t1n<tJv(1) % need to extrapolate to the left (earlier)
                extX1 = t1n:intjsX(end);
                [extjTrj1] = extm(intjsX,extX1,intjTrj,'linear'); % linear extrapolation
                jTrj1 = extjTrj1(:,1:find(extX1==tEr));
            else
                jTrj1 = intjTrj(:,find(intjsX==t1n):find(intjsX==tEr));
            end
            % bin jTrj1
            ss(t).jTrjB = jTrj1(:,1:binSize:size(jTrj1,2)); % 20 ms bins
            % get the joystick set position from jkvt
            ss(t).jsXYbot = jkvt(t).jsTreachPosB(1:2); % joystick bottom
            % compute the reach angle on the horizontal X-Y plane
            ss(t).rchAngDeg = computeReachAngle({ss(t).hTrjB(1:2,:)}, ss(t).jsXYbot);
        end

        %% stim trial info
        if isfield(jkvt,'stimLaserOn')
            if ~isnan(jkvt(t).stimLaserOn) && ~isnan(jkvt(t).stimLaserOff)
                ss(t).tLaserStart = jkvt(t).stimLaserOn;
                ss(t).tLaserStop = jkvt(t).stimLaserOff;
                ss(t).spkTimeBlaserI = ss(t).tLaserStart<=ss(t).spkTimeBins & ss(t).spkTimeBins<=ss(t).tLaserStop;
            end
        end
    end
    fprintf('processed trial %d\n', t) % report trial progression
end

% transfer block information from jkvt to ss
[ss(:).blNumber] = deal(jkvt(:).blNumber);
[ss(:).blType] = deal(jkvt(:).blType);
[ss(:).blShiftLogic] = deal(jkvt(:).blShiftLogic);

% variables to be saved to build trial-type classifiers
ss_trialType = {ss.trialType};
ss_blNumber = {ss.blNumber};
jkvt_rewarded = {jkvt.rewarded};


saveName = filePath(end-19:end-9);

% Save relevant variable for classification analysis
if isfield(ss, 'unitTimeBCtx')
    ss_unitTimeBCtx = {ss.unitTimeBCtx};
    save(fullfile(filePath, strcat('unitTimeBCtx', '_', saveName)), 'ss_unitTimeBCtx')
end

if isfield(ss, 'unitTimeBStr')
    ss_unitTimeBStr = {ss.unitTimeBStr};
    save(fullfile(filePath, strcat('unitTimeBStr', '_', saveName)), 'ss_unitTimeBStr')
end

if isfield(ss, 'unitTimeBCg')
    ss_unitTimeBCg = {ss.unitTimeBCg};
    save(fullfile(filePath, strcat('unitTimeBCg', '_', saveName)), 'ss_unitTimeBCg')
end

if isfield(ss, 'trialType')
    ss_trialType = {ss.trialType};
    save(fullfile(filePath, strcat('trialType', '_', saveName)), 'ss_trialType')
end

if isfield(ss, 'blNumber')
    ss_blNumbs = {ss.blNumber};
    save(fullfile(filePath, strcat('blockNums', '_', saveName)), 'ss_blNumbs')
end

save(fullfile(filePath,strcat('js2p0_tbytSpkHandJsTrjBin_50ms_stimPstim_',saveName)),'ss','jkvt','trI','spkTimesCell*','depth*')

% save(fullfile(filePath,strcat('js2p0_tbytSpkHandJsTrjBin_',saveName)),'depthCtx','depthStr','-append')
% save(fullfile('/Users/parkj/Dropbox (HHMI)/j2p0_dataShare/js2p0_tbytSpkHandJsTrjBin_WR40_081919.mat'),'depthCtx','depthStr','-append')
%load(fullfile(filePath,strcat('js2p0_tbytSpkHandJsTrjBin_',saveName)),'ss','jkvt','trI','spkTimesCell')

%% %%%%%%%%%%%%%%%%%%%
%%% Helper function %%
%%%%%%%%%%%%%%%%%%%%%%
% parse trial/block information
    function jkvt = jkvtBlockParse(jkvt)
        tqd = diff([jkvt(1).pull_torque, jkvt(:).pull_torque]'); % torque change
        p1d = diff([jkvt(1).reachP1, jkvt(:).reachP1]'); % position 1 change

        blNumb = 1;  % block number
        blType = []; % block type
        % detect and parse block shifts across trials
        for tt = 1:size(jkvt,2)
            % first trial
            if tt==1
                jkvt(tt).blNumber = blNumb;
                jkvt(tt).blType = 'first';
                jkvt(tt).blShiftLogic = true;
            end
            % second trial and onward
            if tt>=2
                if tqd(tt)==0 && p1d(tt)==0 % if there's no change
                    jkvt(tt).blNumber = blNumb; % remains the same
                    jkvt(tt).blShiftLogic = false; % not shifted
                    jkvt(tt).blType = []; % block shift type
                else % if there's any change, classify the block shift type
                    blNumb = blNumb+1;
                    jkvt(tt).blNumber = blNumb; % remains the same
                    jkvt(tt).blShiftLogic = true; % shifted
                    % parse torque shift
                    if tqd(tt)>0
                        tqi = 'tqUp';
                    elseif tqd(tt)<0
                        tqi = 'tqDown';
                    elseif tqd(tt)==0
                        tqi = 'tqSame';
                    end
                    % parse position shift
                    if p1d(tt)>0
                        p1i = 'rightward';
                    elseif p1d(tt)<0
                        p1i = 'leftward';
                    elseif p1d(tt)==0
                        if jkvt(tt).reachP1==min([jkvt(:).reachP1])
                            p1i = 'sameleft';
                        elseif jkvt(tt).reachP1==max([jkvt(:).reachP1])
                            p1i = 'sameright';
                        end
                    end
                    jkvt(tt).blType = strcat(tqi,p1i); % block shift type
                end
            end
        end
    end

% get interpolated hand trajectory
    function [intV,xq] = interpsm(x,v)
        xq = x(1):x(end); % new timescale
        inthTrjF = @(a) interp1(x,a,xq); % interpolation function
        vC = mat2cell(v,[ 1 1 1 ], size(v,2)); % convert to cell
        intVC = cellfun(@(a) inthTrjF(a), vC, 'un', 0); % interpolated hTrj cell

        lengthIntVC = unique(cellfun(@length, intVC));

        % to use sg filter get sgfiltFramelen
        if lengthIntVC >= 201
            sgfiltFramelen = 101;
        elseif lengthIntVC < 201
            if mod(lengthIntVC,2)==0
                sgfiltFramelen = lengthIntVC-1; % the frame length for sg filter needs to be an odd number
            else
                sgfiltFramelen = lengthIntVC;
            end
        end
        intVCf = cellfun(@(a) sgolayfilt(a,3,sgfiltFramelen), intVC, 'un', 0); % filtered hTjr cell
        intV = cell2mat(intVCf);
    end

% compute reach angle
    function reachAngle = computeReachAngle( hTrjC, jsXYpos )
        %Angle is calculated relative to the straight line projected from the
        % initial hand position
        distToJs = cellfun(@(a) sum((a-repmat(jsXYpos,1,size(a,2))).^2,1), hTrjC, 'un', 0);  % compute the distance from the jsXY
        for j = 1:length(distToJs)
            [~,tmpMinI] = min(distToJs{j});
            u=hTrjC{j}(:,tmpMinI)-hTrjC{j}(:,1); % reach vector
            %v=[hTrjC{j}(1,tmpMinI);hTrjC{j}(2,1)]-hTrjC{j}(:,1); % reference vector
            v=[hTrjC{j}(1,1);hTrjC{j}(2,tmpMinI)]-hTrjC{j}(:,1); % reference vector
            reachAngle(j) = angleTwoVectors(u,v); % get the angle between reach and reference vectors
        end
    end

    function [targinref] = targetintoreftime(A,timeA,timeRef, varargin)
        %This function takes a times series 'A' across a time window ('timeA') and
        % plug it into the new time window 'timeRef'.
        % 1-ms time resolution is assumed.
        p = targetintoreftime( A,timeA,timeRef, varargin );

        Aval = A(timeRef(1)<=timeA&timeA<=timeRef(end)); % valid portion of A to be placed in timeRef

        if p.Results.fillin == 0
            targinref = zeros(size(A,1),size(timeRef,2));
        elseif isnan(p.Results.fillin)
            targinref = nan(size(A,1),size(timeRef,2));
        end

        tRef1 = find(timeRef>=timeA(1),1,'first'); % 1st point
        if ~isempty(tRef1) && ~isempty(Aval)
            targinref(tRef1:tRef1+min(length(targinref),length(Aval))-1) = Aval;
        end

        function p = targetintoreftime( A,timeA,timeRef, vargs )
            % parse input, and extract name-value pairs
            default_fillin = 0; % by default do not re-read the raw bin file, if done already

            p = inputParser; % create parser object
            addRequired(p,'A');
            addRequired(p,'timeA');
            addRequired(p,'timeRef');
            addParameter(p,'fillin',default_fillin);

            parse(p,A,timeA,timeRef,vargs{:})
        end

    end

end
