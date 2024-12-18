function preprocessingForDecoder_KalmanFilterHTrjF_reachPull(filePath, saveName)
%This function performs preprocessing for state decoding from cortical and
% striatal spike activity using Kalman filter

%% 1. Load data
%clc; clearvars; close all;
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)
% neural data
spkDir = dir('binSpkCountSTRCTX*'); 
load(fullfile(spkDir(1).folder, spkDir(1).name),'spkTimesCell', 'rStartToPull')
%load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','binSpkCountSTRCTXWR40_081919.mat'), 'spkTimesCell', 'rStartToPull', 'p')
S=rStartToPull; clearvars rStartToPull
% behavioral data

behDir = dir('js2p0_tbytSpkHandJsTrjBin_WR*'); 
load(fullfile(behDir(1).folder,fullfile(behDir(1).name)),'jkvt'); 

%% align hand trajectories to neural data
ts = S.currEvt{1}(:,1); % time stamps (e.g. rStartToPull)
% locate trials map each event (ts) to trials in jkvt
tsMapJkvt1 = find(repmat(ts(1),[size(jkvt,2),1]) <= [jkvt(:).trEnd]',1,'first'); 
tsMapJkvt = arrayfun(@(a) find([jkvt(1:end-1).trEnd]'<= repmat(a,[size(jkvt,2)-1,1]) & repmat(a,[size(jkvt,2)-1,1]) <= [jkvt(2:end).trEnd]')+1, ts(2:end), 'un', 0 ); 
tsJkvtTrs = [{tsMapJkvt1}; tsMapJkvt]; 

vfT = {jkvt(:).vFrameTime}'; % video frame time
hTrj = {jkvt(:).hTrjF}'; % hand trajectory
%hTrj = {jkvt(:).hTrjDstBaseXyz}'; % baseline-subtracted hand trajectory

spikeB = S.params.binEdges(2:end);
unitTimeTrial = S.unitTimeTrial;

if ~(length(ts)==size(unitTimeTrial,3))
    error('trial numbers do not match!')
end
spk_kernel = TNC_CreateGaussian(250,15,500,1); % mu, sigma, time, dT

binSize = 20; % 20 ms

% pull torque & reach position combinations
pTqs = round(unique([jkvt(:).pull_torque])./10); % pull torque list (divided by 10 for ordering) 
rPos = round(unique([jkvt(:).reachP1])*10); % reach position list (multiplied by 10 for ordering)
posTq = unique(repmat(pTqs,[length(rPos) 1])+repmat(rPos',[1 length(pTqs)])); % position-torque combinations
posTqCnt = zeros(length(posTq),1); % count occurrence of each combination

% cortex/striatum index
depth = cell2mat(cellfun(@(a) a(2), spkTimesCell(4,:),'un',0))'; % depth from pial surface
ctxI = depth<2000; % <1900; % cortex index
strI = depth>2000; % >2100; % striatum index

s.ctxI = ctxI; 
s.strI = strI; 
s.ctxI_bsc = find(ctxI); 
s.strI_bsc = find(strI); 

% preprocessing of hand trajectories
for t = 1:length(ts)
    % locate the video frames corresponding to current evt
    vfI = cell2mat(cellfun(@(a) a(1)<=ts(t) & a(end)>=ts(t), vfT, 'un', 0)); % vFrame index
    
    if sum(vfI)==1 && ~isempty(hTrj{vfI}) % if hTrj available    
        
        % interpolation of hTrj (time resolution: 4 to 1ms)
        x = vfT{vfI}; % default timescale (4 ms)
        v = hTrj{vfI}; % hTrj to be interpolated
        xq = x(1):x(end); % new timescale
        inthTrj = @(a) interp1(x,a,xq); % interpolation function
        vC = mat2cell(v,[ 1 1 1 ], size(v,2)); % convert to cell
        intVC = cellfun(@(a) inthTrj(a), vC, 'un', 0); % interpolated hTrj cell
        intVCf = cellfun(@(a) sgolayfilt(a,3,201), intVC, 'un', 0); % filtered hTjr cell 
        intV = cell2mat(intVCf); 
        
        % neural spike timepoints
        spikeT = spikeB+ts(t); % timestamps for unitTimeTrial mat (spike counts at 1ms)
        
        t1 = max(spikeT(1), xq(1)); % beh neural common starting point
        tE = min(spikeT(end), xq(end)); % beh neural common end point
        
        t1R = max(t1,ts(t)-100); % 100-ms before rStartToPull
        
        hPos = intV(:,t1<=xq & xq<=tE); % hand position between t1 & tE
        spk  = unitTimeTrial(:,t1<=spikeT & spikeT<=tE ,t); % spk mat between t1 & tE
        
        hPosR = intV(:,t1R<=xq & xq<=tE); % from -100ms relative to reach start
        spkR = unitTimeTrial(:,t1R<=spikeT & spikeT<=tE ,t); % from -100ms relative to reach Start
        
        % force applied onto the joystick 
        tmpForce = jkvt(vfI).movKins.forceMN(1,:); % current trial's smoothed Js force trace (from the encoder, acceleration + Mass)
        tmpForceN = zeros(1, length(tmpForce)); % current trial's force trace only to include negative (means pull) portions 
        tJsNI = jkvt(vfI).movKins.forceMN(1,:)<0; % time index where force applied toward the pull direction
        tmpForceN(tJsNI) = tmpForce(tJsNI); % pull force trace
        sTmpForceN = conv(tmpForceN,spk_kernel,'same'); % smoothing with convolution 
        
        tmpForceNt1RtE = zeros(1,size(hPosR,2)); % force generated onto the joystick in the pull direction 
        tJs = jkvt(vfI).trJsReady:jkvt(vfI).trEnd-1; % time points for the joystick trace from encoder
        tmpT1R = find(tJs==t1R); % t1R on joystick trace
        tmpForceNt1 = sTmpForceN(tmpT1R:min(tmpT1R+size(hPosR,2),length(sTmpForceN))); 
        tmpForceNt1RtE(1:length(tmpForceNt1)) = tmpForceNt1; % length matched force trace (pull portion)
        
        % bundle variables into a structure
        s.spk{1,t} = spk; 
        s.hPos{1,t} = hPos; % hand position xyz 
        s.hVel{1,t} = diff([hPos(:,1) hPos],1,2)./(1/1000)./10; % cm/s
        
        s.hPosR{1,t} = hPosR; % hand position xyz cut from -100ms relative to reach Start 
        s.spkR{1,t} = spkR; % spike count mat cut from -100ms relative to reach Start 
        s.jsForceR{1,t} = tmpForceNt1RtE; % force trace cut from -100ms relative to reach Start
        
        % timestamps all in ms
        s.time(t).t1 = t1; % start time point of spk and hPos
        s.time(t).tE = tE; % end time point of spk and hPos
        s.time(t).t1R = t1R; % start time point of spkR and hPosR (cut from -100-ms relative to reachStart)
        s.time(t).ts = ts(t); % reference time point of the current trial (e.g. rStartToPull)
        s.time(t).jkvtTr = tsJkvtTrs{t}; % trial # in jkvt
             
        %% bin spkR and hPosR to construct dat
        % identify pull torque and reach position
        posTqI = find(jkvt(tsJkvtTrs{t}).pull_torque/10 + jkvt(tsJkvtTrs{t}).reachP1*10 == posTq); 
        posTqCnt(posTqI)=posTqCnt(posTqI)+1; 
        
        r = posTqCnt(posTqI); % row index
        c = posTqI; % column index
        
        s.dat.pos1{r,c} = jkvt(tsJkvtTrs{t}).reachP1; % store reach position
        s.dat.trq{r,c} = jkvt(tsJkvtTrs{t}).pull_torque; % store pull torque
        
        timeBin = binSize:binSize:size(hPosR,2); 
        curhPosR = hPosR(:,1:timeBin(end)); 
        curBpos = curhPosR(:,timeBin); 
        curBvel = (curBpos(:,2:end)-curBpos(:,1:end-1))*(1000/binSize)/10; % velocity (cm/s) 
        
        curForceR = tmpForceNt1RtE(:,1:timeBin(end));
        curBforce = curForceR(timeBin); 
        
        curSpkR = spkR(:,1:timeBin(end-1))'; % timebin(1ms)-by-neuron spike matrix 
        rsCurSpkR = reshape(curSpkR,binSize,[],size(curSpkR,2)); 
        binSpkRCtxStr = squeeze(sum(rsCurSpkR,1))'; % unit-by-timeBin        
        
        % get the pullStart and pullStop point and pull index, if exists
        if isfield(jkvt(vfI).movKins,'pullStart') && isfield(jkvt(vfI).movKins,'pullStop')
            pullStart = jkvt(vfI).trJsReady + jkvt(vfI).movKins.pullStart;
            pullStop = jkvt(vfI).trJsReady + jkvt(vfI).movKins.pullStop;
            
            s.time(t).tPullStart = pullStart; 
            s.time(t).tPullStop = pullStop; 
            
            t1tEpullI = pullStart<=t1:tE & t1:tE<=pullStop;
            t1RtEpullI = pullStart<=t1R:tE & t1R:tE<=pullStop;
            
            s.dat.pullIdx{r,c} = sum(reshape(t1RtEpullI(1:timeBin(end-1)), binSize, []))>=1; 
            tmpPullStart = max(1,find(s.dat.pullIdx{r,c},1,'first')); 
            tmpPullEnd = max(tmpPullStart,find(s.dat.pullIdx{r,c},1,'last')); 
            
            tmpState = [curBpos(:,1:end-1); curBvel(:,1:end); curBforce(:,1:end-1)]; % state 3-d pos; 3-d vel  
            
            if ~isempty(tmpPullStart) && ~isempty(tmpPullEnd)
                s.dat.stateR{r,c} = tmpState(:,1:tmpPullStart-1);  
                s.dat.spkCtxR{r,c} = binSpkRCtxStr(ctxI,1:tmpPullStart-1); 
                s.dat.spkStrR{r,c} = binSpkRCtxStr(strI,1:tmpPullStart-1);      
                s.dat.spkCtxStrR{r,c} = binSpkRCtxStr(:,1:tmpPullStart-1);            
            
                s.dat.stateP{r,c} = tmpState(:,tmpPullStart:min(size(tmpState,2),tmpPullEnd+5));  
                s.dat.spkCtxP{r,c} = binSpkRCtxStr(ctxI,tmpPullStart:min(size(tmpState,2),tmpPullEnd+5)); 
                s.dat.spkStrP{r,c} = binSpkRCtxStr(strI,tmpPullStart:min(size(tmpState,2),tmpPullEnd+5));      
                s.dat.spkCtxStrP{r,c} = binSpkRCtxStr(:,tmpPullStart:min(size(tmpState,2),tmpPullEnd+5));    
            end
            s.dat.trialEvt{r,c} = t;
            s.dat.trialJkvt{r,c} = tsJkvtTrs{t};            
        end
        % get the laserStart and laserStop point and laser index, if exists
        if isfield(jkvt,'stimLaserOn')
            if ~isempty(jkvt(vfI).stimLaserOn) && ~isempty(jkvt(vfI).stimLaserOff)
                s.time(t).tLaserStart = jkvt(vfI).stimLaserOn;
                s.time(t).tLaserStop = jkvt(vfI).stimLaserOff;
                
                t1tElaserI = s.time(t).tLaserStart<=t1:tE & t1:tE<=s.time(t).tLaserStop;
                t1RtElaserI = s.time(t).tLaserStart<=t1R:tE & t1R:tE<=s.time(t).tLaserStop;
                
                s.dat.laserOffTime{r,c} = s.time(t).tLaserStop-t1R; % when did laser went off relative to the reach start
                
                if isfield(jkvt(vfI).movKins,'pullStart') && isfield(jkvt(vfI).movKins,'pullStop')
                      s.dat.laserIdx{r,c} = sum(reshape(t1RtElaserI(1:timeBin(end-1)), binSize, []))>=1;
                      s.dat.laserIdxR{r,c} = s.dat.laserIdx{r,c}(:,1:tmpPullStart-1); 
                      s.dat.laserIdxP{r,c} = s.dat.laserIdx{r,c}(:,tmpPullStart:min(size(tmpState,2),tmpPullEnd+5)); 
                end
            end
        else
            s.dat.laserIdx{r,c} = 0; 
            s.dat.laserIdxR{r,c} = 0; 
            s.dat.laserIdxP{r,c} = 0; 
        end
    end
    fprintf('processed event# %d\n', t)  
end
save(fullfile(filePath,strcat('preprocessKFdecodeHTrjCtxStr_reachPull_',saveName)),'s')
end