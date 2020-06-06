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

% preprocessing of hand trajectories
for t = 1:length(ts)
    % locate the video frames corresponding to current evt
    vfI = cell2mat(cellfun(@(a) a(1)<=ts(t) & a(end)>=ts(t), vfT, 'un', 0)); % vFrame index
    
    if sum(vfI)==1 && ~isempty(hTrj{vfI}) % if hTrj available    
        % interpolation
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
        
        hPos = intV(:,t1<=xq & xq<=tE);
        spk  = unitTimeTrial(:,t1<=spikeT & spikeT<=tE ,t);
        
        % get the pullStop point, if exists
        if isfield(jkvt(vfI).movKins,'pullStart') && isfield(jkvt(vfI).movKins,'pullStop')
            pullStart = jkvt(vfI).trJsReady + jkvt(vfI).movKins.pullStart;
            pullStop = jkvt(vfI).trJsReady + jkvt(vfI).movKins.pullStop;
            
            s.time(t).tPullStart = pullStart; 
            s.time(t).tPullStop = pullStop; 
        end
        
        % build input and target structures
        conv_spk = zeros(size(spk,1), size(spk,2)); 
        for u = 1:size(spk,1)
            conv_spk(u,:) = conv(spk(u,:), spk_kernel, 'same');
        end
        
        s.spk{1,t} = conv_spk; 
        s.hPos{1,t} = hPos; % hand position xyz 
        s.hVel{1,t} = diff([hPos(:,1) hPos],1,2)./(1/1000)./10; % cm/s
        
        s.time(t).t1 = t1; 
        s.time(t).tE = tE; 
        
        s.ts(t,1) = ts(t); % time in ms
        s.ts(t,2) = cell2mat(tsJkvtTrs(t)); % trial # in jkvt
    end
end

s.spk_strIdx = cell2mat(S.isStr); 
save(fullfile(filePath,'preprocessDecodeHTrjCtxStr_WR40_081919.mat'),'s')
