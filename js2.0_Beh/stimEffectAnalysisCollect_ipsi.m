% pullStart latency distribution
filePathIT = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR46_021420'};
filePathPT = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR43_021620'};


%% IT
for f = 1:length(filePathIT)
    cd(filePathIT{f})
    
    if exist(fullfile(filePathIT{f},'jsTime1k_Kinematics_VideoFiles.mat'),'file')
        S = load(fullfile(filePathIT{f},'jsTime1k_Kinematics_VideoFiles.mat'));  % load jsTime1k_KV.mat
        S = S.('jsTime1k_KV');
    elseif exist(fullfile(filePathIT{f},'jsTime1k_Kinematics.mat'),'file')
        S = load(fullfile(filePathIT{f},'jsTime1k_Kinematics.mat'));  % load jsTime1k_KV.mat
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
    
    for t = 1:length(S)
        if S(t).rewarded && strcmpi(S(t).trialType,'sp')
            if isnan(S(t).stimLaserOn)
                if isstruct(S(t).movKins)
                    if isfield(S,'pLaserOn') && ~isnan(S(t).pLaserOn) % if pseudo laser pulse was used
                        trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - S(t).pLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
                    end
                end
            elseif ~isnan(S(t).stimLaserOn)
                if isstruct(S(t).movKins)
                    trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - S(t).stimLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
                end
            end
        end
    end
    clearvars t
    
    rStartRelITStim.stmhTq{f,1} = trStartPt(stmTrials & hTqTrs);
    rStartRelITStim.stmlTq{f,1} = trStartPt(stmTrials & lTqTrs);
    rStartRelITStim.pStmhTq{f,1} = trStartPt(~stmTrials & hTqTrs);
    rStartRelITStim.pStmlTq{f,1} = trStartPt(~stmTrials & lTqTrs);
   
    fprintf('processed file %d\n', f) % report unit progression
end

% IT high torque
rStart.stimIThTq = cell2mat(rStartRelITStim.stmhTq(:,1));
rStart.stimIThTq = rStart.stimIThTq(rStart.stimIThTq<10000 & rStart.stimIThTq>0);

rStart.pStimIThTq = cell2mat(rStartRelITStim.pStmhTq(:,1)); 
rStart.pStimIThTq = rStart.pStimIThTq(rStart.pStimIThTq<10000 & rStart.pStimIThTq>0); 

% IT low torque
rStart.stimITlTq = cell2mat(rStartRelITStim.stmlTq(:,1)); 
rStart.stimITlTq = rStart.stimITlTq(rStart.stimITlTq<10000 & rStart.stimITlTq>0);

rStart.pStimITlTq = cell2mat(rStartRelITStim.pStmlTq(:,1)); 
rStart.pStimITlTq = rStart.pStimITlTq(rStart.pStimITlTq<10000 & rStart.pStimITlTq>0); 

%% histogram the pullStart latencies 
edges = 0:250:10000; 
% IT low torque pStim
figure; hold on; 
rStart.pStimITlTqHc = histcounts(rStart.pStimITlTq, edges);
rStart.pStimITlTqHcNorm = rStart.pStimITlTqHc./sum(rStart.pStimITlTqHc);
bar(edges(1:end-1),rStart.pStimITlTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

% IT low torque stim
rStart.stimITlTqHc = histcounts(rStart.stimITlTq, edges);
rStart.stimITlTqHcNorm = rStart.stimITlTqHc./sum(rStart.stimITlTqHc);
bar(edges(1:end-1),rStart.stimITlTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
yticks(0:.1:1);
set(gca,'TickDir','out')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','pullStartTime_ITlowTq_pStim_Stim_ipsi'),'-dpdf','-painters')

% IT high torque pStim
figure; hold on; 
rStart.pStimIThTqHc = histcounts(rStart.pStimIThTq, edges);
rStart.pStimIThTqHcNorm = rStart.pStimIThTqHc./sum(rStart.pStimIThTqHc);
bar(edges(1:end-1),rStart.pStimIThTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

% IT high torque stim
rStart.stimIThTqHc = histcounts(rStart.stimIThTq, edges);
rStart.stimIThTqHcNorm = rStart.stimIThTqHc./sum(rStart.stimIThTqHc);
bar(edges(1:end-1),rStart.stimIThTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 .6]);
yticks(0:.1:1);
set(gca,'TickDir','out')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','pullStartTime_IThighTq_pStim_Stim_ipsi'),'-dpdf','-painters')

%% PT
for f = 1:length(filePathPT)
    cd(filePathPT{f})
    
    if exist(fullfile(filePathPT{f},'jsTime1k_Kinematics_VideoFiles.mat'),'file')
        S = load(fullfile(filePathPT{f},'jsTime1k_Kinematics_VideoFiles.mat'));  % load jsTime1k_KV.mat
        S = S.('jsTime1k_KV');
    elseif exist(fullfile(filePathPT{f},'jsTime1k_Kinematics.mat'),'file')
        S = load(fullfile(filePathPT{f},'jsTime1k_Kinematics.mat'));  % load jsTime1k_KV.mat
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
    
    for t = 1:length(S)
        if S(t).rewarded && strcmpi(S(t).trialType,'sp')
            if isnan(S(t).stimLaserOn)
                if isstruct(S(t).movKins)
                    if isfield(S,'pLaserOn') && ~isnan(S(t).pLaserOn) % if pseudo laser pulse was used
                        trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - S(t).pLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
                    end
                end
            elseif ~isnan(S(t).stimLaserOn)
                if isstruct(S(t).movKins)
                    trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - S(t).stimLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
                end
            end
        end
    end
    clearvars t
    
    rStartRelPTStim.stmhTq{f,1} = trStartPt(stmTrials & hTqTrs);
    rStartRelPTStim.stmlTq{f,1} = trStartPt(stmTrials & lTqTrs);
    rStartRelPTStim.pStmhTq{f,1} = trStartPt(~stmTrials & hTqTrs);
    rStartRelPTStim.pStmlTq{f,1} = trStartPt(~stmTrials & lTqTrs);
   
    fprintf('processed file %d\n', f) % report unit progression
end

% PT high torque
rStart.stimPThTq = cell2mat(rStartRelPTStim.stmhTq(:,1));
rStart.stimPThTq = rStart.stimPThTq(rStart.stimPThTq<10000 & rStart.stimPThTq>0);

rStart.pStimPThTq = cell2mat(rStartRelPTStim.pStmhTq(:,1)); 
rStart.pStimPThTq = rStart.pStimPThTq(rStart.pStimPThTq<10000 & rStart.pStimPThTq>0); 

% PT low torque
rStart.stimPTlTq = cell2mat(rStartRelPTStim.stmlTq(:,1)); 
rStart.stimPTlTq = rStart.stimPTlTq(rStart.stimPTlTq<10000 & rStart.stimPTlTq>0);

rStart.pStimPTlTq = cell2mat(rStartRelPTStim.pStmlTq(:,1)); 
rStart.pStimPTlTq = rStart.pStimPTlTq(rStart.pStimPTlTq<10000 & rStart.pStimPTlTq>0); 

%% histogram the pullStart latencies 
edges = 0:250:10000; 
% PT low torque pStim
figure; hold on; 
rStart.pStimPTlTqHc = histcounts(rStart.pStimPTlTq, edges);
rStart.pStimPTlTqHcNorm = rStart.pStimPTlTqHc./sum(rStart.pStimPTlTqHc);
bar(edges(1:end-1),rStart.pStimPTlTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

% PT low torque stim
rStart.stimPTlTqHc = histcounts(rStart.stimPTlTq, edges);
rStart.stimPTlTqHcNorm = rStart.stimPTlTqHc./sum(rStart.stimPTlTqHc);
bar(edges(1:end-1),rStart.stimPTlTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 .4])
yticks(0:.1:1);
set(gca,'TickDir','out')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','pullStartTime_PTlowTq_pStim_Stim_ipsi'),'-dpdf','-painters')

% PT high torque pStim
figure; hold on; 
rStart.pStimPThTqHc = histcounts(rStart.pStimPThTq, edges);
rStart.pStimPThTqHcNorm = rStart.pStimPThTqHc./sum(rStart.pStimPThTqHc);
bar(edges(1:end-1),rStart.pStimPThTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

% PT high torque stim
rStart.stimPThTqHc = histcounts(rStart.stimPThTq, edges);
rStart.stimPThTqHcNorm = rStart.stimPThTqHc./sum(rStart.stimPThTqHc);
bar(edges(1:end-1),rStart.stimPThTqHcNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 .4])
yticks(0:.1:1);
set(gca,'TickDir','out')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','pullStartTime_PThighTq_pStim_Stim_ipsi'),'-dpdf','-painters')












