%This script examines the behavioral effect of contra-M1 silencing with 1) latency for pullStart  

% pullStart latency distribution
filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_090319_m1',...
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_050119'...
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081519'};

%% 
for f = 1:length(filePath)
    cd(filePath{f})
    
    if exist(fullfile(filePath{f},'jsTime1k_Kinematics_VideoFiles.mat'),'file')
        S = load(fullfile(filePath{f},'jsTime1k_Kinematics_VideoFiles.mat'));  % load jsTime1k_KV.mat
        S = S.('jsTime1k_KV');
    elseif exist(fullfile(filePath{f},'jsTime1k_Kinematics.mat'),'file')
        S = load(fullfile(filePath{f},'jsTime1k_Kinematics.mat'));  % load jsTime1k_KV.mat
        S = S.('jsTime1k_K');
    else
        error('No behavioral data found!!!')
    end
     
    tqs = unique([S(:).pull_torque]); % pull torques in the session
    hTqTrs = cellfun(@ (a) a==tqs(end), {S(:).pull_torque})'; % high torque trials
    lTqTrs = cellfun(@ (a) a==tqs(1), {S(:).pull_torque})'; % low torque trials
    stmTrials = ~isnan([S(:).stimLaserOn]');
    
    lTqPullStart = [];
    hTqPullStart = [];
    
    trStartPt = nan(length(S),1);
    trJsReadyRelpLaser = nan(length(S),1); 
    trJsReadyRelLaser = nan(length(S),1); 
    
    for t = 1:length(S)
        if S(t).rewarded && strcmpi(S(t).trialType,'sp')
            if isnan(S(t).stimLaserOn)
                if isstruct(S(t).movKins)
                    if isfield(S,'pLaserOn') && ~isnan(S(t).pLaserOn) % if pseudo laser pulse was used
                        trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - (S(t).pLaserOn+4000); % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
                        trJsReadyRelpLaser(t,1) = S(t).trJsReady-S(t).pLaserOn; % time difference
                    end
                end
            elseif ~isnan(S(t).stimLaserOn)
                if isstruct(S(t).movKins)
                    trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - (S(t).stimLaserOn+4000); % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
                    trJsReadyRelLaser(t,1) = S(t).trJsReady-S(t).stimLaserOn; % time difference
                end
            end
        end
    end
    clearvars t
    
    rStartRelStim.stmhTq{f,1} = trStartPt(stmTrials & hTqTrs);
    rStartRelStim.stmlTq{f,1} = trStartPt(stmTrials & lTqTrs);
    rStartRelStim.pStmhTq{f,1} = trStartPt(~stmTrials & hTqTrs);
    rStartRelStim.pStmlTq{f,1} = trStartPt(~stmTrials & lTqTrs);
    rStartRelStim.trJsReadyRelpLaser{f,1} = trJsReadyRelpLaser; 
    rStartRelStim.trJsReadyRelLaser{f,1} = trJsReadyRelLaser; 
    
    fprintf('processed file %d\n', f) % report unit progression
end

% high torque
rStart.stimhTq = cell2mat(rStartRelStim.stmhTq(:,1));
rStart.stimhTq = rStart.stimhTq(rStart.stimhTq<10000); % & rStart.stimhTq>0);

rStart.pStimhTq = cell2mat(rStartRelStim.pStmhTq(:,1)); 
rStart.pStimhTq = rStart.pStimhTq(rStart.pStimhTq<10000); % & rStart.pStimhTq>0); 

% low torque
rStart.stimlTq = cell2mat(rStartRelStim.stmlTq(:,1)); 
rStart.stimlTq = rStart.stimlTq(rStart.stimlTq<10000); % & rStart.stimlTq>0);

rStart.pStimlTq = cell2mat(rStartRelStim.pStmlTq(:,1)); 
rStart.pStimlTq = rStart.pStimlTq(rStart.pStimlTq<10000); % & rStart.pStimlTq>0); 

%% histogram the pullStart latencies 
edges = -2000:250:10000; 
%  low torque pStim
figure; hold on; 
rStart.pStimlTqHc = histcounts(rStart.pStimlTq, edges);
rStart.pStimlTqHcNorm = rStart.pStimlTqHc./sum(rStart.pStimlTqHc); % normalize
rStart.pStimlTqHcNormCumSum = cumsum(rStart.pStimlTqHcNorm); % cumulative sum
bar(edges(1:end-1),rStart.pStimlTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

%  low torque stim
rStart.stimlTqHc = histcounts(rStart.stimlTq, edges);
rStart.stimlTqHcNorm = rStart.stimlTqHc./sum(rStart.stimlTqHc);
rStart.stimlTqHcNormCumSum = cumsum(rStart.stimlTqHcNorm); % cumulative sum
bar(edges(1:end-1),rStart.stimlTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 .3]);
yticks(0:.1:1);
set(gca,'TickDir','out')
%print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','pullStartTime_lowTq_pStim_Stim_ipsi'),'-dpdf','-painters')

%  high torque pStim
figure; hold on; 
rStart.pStimhTqHc = histcounts(rStart.pStimhTq, edges);
rStart.pStimhTqHcNorm = rStart.pStimhTqHc./sum(rStart.pStimhTqHc);
rStart.pStimhTqHcNormCumSum = cumsum(rStart.pStimhTqHcNorm); % cumulative sum
bar(edges(1:end-1),rStart.pStimhTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

%  high torque stim
rStart.stimhTqHc = histcounts(rStart.stimhTq, edges);
rStart.stimhTqHcNorm = rStart.stimhTqHc./sum(rStart.stimhTqHc);
rStart.stimhTqHcNormCumSum = cumsum(rStart.stimhTqHcNorm); % cumulative sum
bar(edges(1:end-1),rStart.stimhTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 .3]);
yticks(0:.1:1);
set(gca,'TickDir','out')
%print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','pullStartTime_highTq_pStim_Stim_ipsi'),'-dpdf','-painters')

% plot cumulative pullStart
figure; hold on; 
plot(edges(1:end-1), rStart.pStimhTqHcNormCumSum)
plot(edges(1:end-1), rStart.stimhTqHcNormCumSum)
