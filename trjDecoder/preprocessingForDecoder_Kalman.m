%% 1. Load data
%clc; clearvars; close all;

filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles';
cd(filePath)
% neural data
load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','binSpkCountSTRCTXWR40_081919.mat'), 'spkTimesCell', 'rStartToPull', 'p')
S=rStartToPull; clearvars rStartToPull
% behavioral data
load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','jsTime1k_KinematicsTrajectories.mat'),'jkvt','meta')

%% align hand trajectories to neural data
ts = S.currEvt{1}(:,1); % time stamps (e.g. rStartToPull)
% locate trials map each event (ts) to trials in jkvt
tsMapJkvt1 = find(repmat(ts(1),[size(jkvt,2),1]) <= [jkvt(:).trEnd]',1,'first'); 
tsMapJkvt = arrayfun(@(a) find([jkvt(1:end-1).trEnd]'<= repmat(a,[size(jkvt,2)-1,1]) & repmat(a,[size(jkvt,2)-1,1]) <= [jkvt(2:end).trEnd]')+1, ts(2:end), 'un', 0 ); 
tsJkvtTrs = [{tsMapJkvt1}; tsMapJkvt]; 

vfT = {jkvt(:).vFrameTime}'; % video frame time
hTrj = {jkvt(:).hTrjDstBaseXyz}'; % baseline-subtracted hand trajectory

spikeB = S.params.binEdges(2:end);
unitTimeTrial = S.unitTimeTrial;

if ~(length(ts)==size(unitTimeTrial,3))
    error('trial numbers do not match!')
end
spk_kernel = TNC_CreateGaussian(500,25,1000,1); % mu, sigma, time, dT

binSize = 20; % 20 ms

% pull torque & reach position combinations
pTqs = round(unique([jkvt(:).pull_torque])./10); % pull torque list (divided by 10 for ordering) 
rPos = round(unique([jkvt(:).reachP1])*10); % reach position list (multiplied by 10 for ordering)
posTq = unique(repmat(pTqs,[length(rPos) 1])+repmat(rPos',[1 length(pTqs)])); % position-torque combinations
posTqCnt = zeros(length(posTq),1); % count occurrence of each combination

% cortex/striatum index
depth = cell2mat(cellfun(@(a) a(2), spkTimesCell(4,:),'un',0))'; % depth from pial surface
ctxI = depth<1900; % cortex index
strI = depth>2100; % striatum index

% preprocessing of hand trajectories
for t = 1:length(ts)
    % locate the video frames corresponding to current evt
    vfI = cell2mat(cellfun(@(a) a(1)<=ts(t) & a(end)>=ts(t), vfT, 'un', 0)); % vFrame index
    
    if sum(vfI)==1 && ~isempty(hTrj{vfI}) % if hTrj available    
        
        % interpolation of hTrj (time resolution: 4 to 1ms)
        x = vfT{vfI}; % default timescale (4 ms)
        v = hTrj{vfI}; % hTrj to be interpolated
        xq = vfT{vfI}(1):vfT{vfI}(end); % new timescale
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
        
        hPos = intV(:,t1<=xq & xq<=tE);
        spk  = unitTimeTrial(:,t1<=spikeT & spikeT<=tE ,t);
        
        hPosR = intV(:,t1R<=xq & xq<=tE); % from -100ms relative to reach start
        spkR = unitTimeTrial(:,t1R<=spikeT & spikeT<=tE ,t); % from -100ms relative to reach Start
        
        s.spk{1,t} = spk; 
        s.hPos{1,t} = hPos; % hand position xyz 
        s.hVel{1,t} = diff([hPos(:,1) hPos],1,2)./(1/1000)./10; % cm/s
        
        s.hPosR{1,t} = hPosR; % hand position xyz cut from -100ms relative to reach Start 
        s.spkR{1,t} = spkR; % spike count mat cut from -100ms relative to reach Start 
        
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
        
        timeBin = binSize:binSize:size(hPosR,2); 
        curhPosR = hPosR(:,1:timeBin(end)); 
        curBpos = curhPosR(:,timeBin); 
        curBvel = (curBpos(:,2:end)-curBpos(:,1:end-1))*(1000/binSize)/10; % velocity (cm/s) 
        s.dat.state{r,c} = [curBpos(:,1:end-1); curBvel]; % state 3-d pos; 3-d vel  
        
        curSpkR = spkR(:,1:timeBin(end-1))'; % timebin(1ms)-by-neuron spike matrix 
        rsCurSpkR = reshape(curSpkR,binSize,[],size(curSpkR,2)); 
        binSpkRCtxStr = squeeze(sum(rsCurSpkR,1))'; % unit-by-timeBin
        
        s.dat.spkCtx{r,c} = binSpkRCtxStr(ctxI,:); 
        s.dat.spkStr{r,c} = binSpkRCtxStr(strI,:);      
        s.dat.spkCtxStr{r,c} = binSpkRCtxStr;            
        s.dat.trialEvt{r,c} = t;
        s.dat.trialJkvt{r,c} = tsJkvtTrs{t}; 
        
        pulldetectI = 0; 
        % get the pullStart and pullStop point, if exists
        if isfield(jkvt(vfI).movKins,'pullStart') && isfield(jkvt(vfI).movKins,'pullStop')
            pulldetectI = 1; 
            pullStart = jkvt(vfI).trJsReady + jkvt(vfI).movKins.pullStart;
            pullStop = jkvt(vfI).trJsReady + jkvt(vfI).movKins.pullStop;
            
            s.time(t).tPullStart = pullStart; 
            s.time(t).tPullStop = pullStop; 
            
            t1tEpullI = pullStart<=t1:tE & t1:tE<=pullStop;
            t1RtEpullI = pullStart<=t1R:tE & t1R:tE<=pullStop;
            
            s.dat.pullIdx{r,c} = sum(reshape(t1RtEpullI(1:timeBin(end-1)), binSize, []))>=1; 
        end
      
    end
    fprintf('processed event# %d\n', t)  
end

save(fullfile(filePath,'preprocessKFdecodeHTrjCtxStr_WR40_081919.mat'),'s')
