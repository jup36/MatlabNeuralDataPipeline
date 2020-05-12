function stimEffectAnalysisJsTrj(filePath)

%filePath = '/Volumes/RAID2/parkj/NeuralData/js2.0/WR37/022619/Matfiles';
cd(filePath)

if exist(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'),'file')
    S = load(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'));  % load jsTime1k_KV.mat
    S = S.('jsTime1k_KV');
elseif exist(fullfile(filePath,'jsTime1k_Kinematics.mat'),'file')
    S = load(fullfile(filePath,'jsTime1k_Kinematics.mat'));  % load jsTime1k_KV.mat
    S = S.('jsTime1k_K');
else
    error('No behavioral data found!!!')
end
    
maxPtsToPlot = 20000; % 10s
rewardIdx = [S.rewarded];
trJsReady = [S.trJsReady];
firstRwdTrial = find([S(:).rewarded]==1,1,'first'); 

%% plot sgJsTrajmm of stim vs non-stim trials aligned to the actual and putative stim onset times
% All trials
figure; hold on;
for t = 2:length(S)
    
    if isnan(S(t).stimLaserOn)
        if isstruct(S(t).movKins)         
            if isfield(S,'pLaserOn') % if pseudo laser pulse was used 
              trStartPt = S(t).trJsReady-S(t).pLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            else
              trStartPt = S(t).trJsReady-(S(t-1).trEnd+4000); % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            end        
            plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[70,240,240]./255);
            plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [70,240,240]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end
    elseif ~isnan(S(t).stimLaserOn)
        if isstruct(S(t).movKins)
            trStartPt = S(t).trJsReady-S(t).stimLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[240,50,230]./255);
            plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [240,50,230]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end     
    end
    
end
clearvars t trStartPt
hold off;
xlim([0 8000]); xlabel('Time relative to actual/putative stim onset')
ylim([-20 20]); ylabel('Joystick position (mm)')
print(fullfile(filePath,'Figure','JsTrajStimVsNoStim_AllTrials'), '-dpdf','-painters')

% Successful trials only
hold on; 
for t = 2:length(S)
    
    if isnan(S(t).stimLaserOn)
        if isstruct(S(t).movKins) && strcmpi(S(t).trialType,'sp')
            if isfield(S,'pLaserOn') % if pseudo laser pulse was used 
              trStartPt = S(t).trJsReady-S(t).pLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
              trStartPtCollect(t,1) = trStartPt; 
            else
              trStartPt = S(t).trJsReady-(S(t-1).trEnd+4000); % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            end   
            plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[70,240,240]./255);
            plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [70,240,240]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end
    elseif ~isnan(S(t).stimLaserOn) && strcmpi(S(t).trialType,'sp')
        if isstruct(S(t).movKins)
            trStartPt = S(t).trJsReady-S(t).stimLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[240,50,230]./255);
            plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [240,50,230]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end     
    end
   
end
clearvars t trStartPt
xlim([0 8000]); xlabel('Time relative to actual/putative stim onset(ms)')
ylim([-20 5]); ylabel('Joystick position (mm)')
print(fullfile(filePath,'Figure','JsTrajStimVsNoStim-SuccessTrials'), '-dpdf','-painters')
hold off;

%% Reach probability P(reachOn|Time), Expected Js position E(JsPosition|Time)
valTrCnt = 0; 
for t = 1:length(S)
  
    if isstruct(S(t).movKins) && t>firstRwdTrial %%&& strcmpi(S(t).trialType,'sp')
        valTrCnt = valTrCnt + 1; 
        tempTrj = S(t).movKins.sgJsTrajmm; 
        % mark reach On before each time point
        jsReachOn{valTrCnt} = zeros(1,10000);
        if strcmpi(S(t).trialType,'sp')
            jsReachOn{valTrCnt}(S(t).movKins.pullStart:end) = 1;
        end
        
        % grab the traj
        if length(tempTrj)<=10000
            jsPosE{valTrCnt}(1:length(tempTrj)) = tempTrj; 
            jsPosE{valTrCnt}(length(tempTrj)+1:10000) = tempTrj(end);             
            
        else 
            jsPosE{valTrCnt} = tempTrj(1:10000);     
        end
            
        % get stim idx
        if isnan(S(t).stimLaserOn)
            stimTrIdx{valTrCnt} = 'noStim';     
        elseif ~isnan(S(t).stimLaserOn)    
            stimTrIdx{valTrCnt} = 'stim'; 
        end
    end

end

X = 1:10000; 

clear g
g(1,1)=gramm('x',X,'y',jsPosE,'color',stimTrIdx);
g(1,2)=gramm('x',X,'y',jsReachOn,'color',stimTrIdx);

g(1,1).stat_summary();
g(1,1).set_title('E(jsPos)');
g(1,1).axe_property('ylim',[-15 1]); 
g(1,1).axe_property('xlim',[0 5000]); 

g(1,2).stat_summary();
g(1,2).set_title('P(reachOn|Time)');
g(1,2).axe_property('xlim',[0 5000]); 

figure('Position',[100 100 800 550]);
g.draw();
%set(gcf,'renderer','painters')
print(fullfile(filePath,'Figure','reachProb_expectedJsPos'), '-dpdf','-bestfit','-painters')

[~,fileN,~] = fileparts(filePath); 
save(fullfile(filePath,strcat(fileN,'_stimR')))

%% plot sgJsTrajmm of stim vs non-stim trials aligned to the actual and putative stim onset times
% All trials
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
                    trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - (S(t).pLaserOn+4000); % align to the putative pStim on, stim generated 4-sec after the last trial offset with the camTrigger                    
                    
                end
            end
        elseif ~isnan(S(t).stimLaserOn)
            if isstruct(S(t).movKins)
                trStartPt(t,1) = S(t).trJsReady + S(t).movKins.pullStart - (S(t).stimLaserOn+4000); % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            end
        end
    end
end
clearvars t

rStartRelStim.stmhTq = trStartPt(stmTrials & hTqTrs); 
rStartRelStim.stmlTq = trStartPt(stmTrials & lTqTrs); 
rStartRelStim.pStmhTq = trStartPt(~stmTrials & hTqTrs); 
rStartRelStim.pStmlTq = trStartPt(~stmTrials & lTqTrs); 

rStartRelStim.stmhTq = rStartRelStim.stmhTq(rStartRelStim.stmhTq<10000);% & rStartRelStim.stmhTq>0); % stim high torque
rStartRelStim.stmlTq = rStartRelStim.stmlTq(rStartRelStim.stmlTq<10000);% & rStartRelStim.stmlTq>0); % stim low torque
rStartRelStim.pStmhTq = rStartRelStim.pStmhTq(rStartRelStim.pStmhTq<10000);% & rStartRelStim.pStmhTq>0); % pStim high torque
rStartRelStim.pStmlTq = rStartRelStim.pStmlTq(rStartRelStim.pStmlTq<10000);% & rStartRelStim.pStmlTq>0); % pStim low torque

%% histogram the pullStart latencies 
edges = -4000:200:10000; 
% low torque stim
figure; hold on; 
[rStartRelStim.stmlTqHist, rStartRelStim.stmlTqHistE] = histcounts(rStartRelStim.stmlTq, edges); % histcounts
rStartRelStim.stmlTqHistNorm = rStartRelStim.stmlTqHist./sum(rStartRelStim.stmlTqHist); % normalized histcounts 
rStartRelStim.stmlTqHcNormCumSum = cumsum(rStartRelStim.stmlTqHistNorm); % cumulative sum
bar(edges(1:end-1),rStartRelStim.stmlTqHistNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 

% low torque pStim
[rStartRelStim.pStmlTqHist, rStartRelStim.pStmlTqHistE] = histcounts(rStartRelStim.pStmlTq, edges);
rStartRelStim.pStmlTqHistNorm = rStartRelStim.pStmlTqHist./sum(rStartRelStim.pStmlTqHist); 
rStartRelStim.pStmlTqHcNormCumSum = cumsum(rStartRelStim.pStmlTqHistNorm); % cumulative sum
bar(edges(1:end-1),rStartRelStim.pStmlTqHistNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 .4])
yticks(0:.1:1);
set(gca,'TickDir','out')

% plot cumulative pullStart on low torque
figure; hold on; 
plot(edges(1:end-1), rStartRelStim.pStmlTqHcNormCumSum)
plot(edges(1:end-1), rStartRelStim.stmlTqHcNormCumSum)
set(gca,'TickDir','out')

% high torque stim
figure; hold on; 
[rStartRelStim.stmhTqHist, rStartRelStim.stmhTqHistE] = histcounts(rStartRelStim.stmhTq, edges);
rStartRelStim.stmhTqHistNorm = rStartRelStim.stmhTqHist./sum(rStartRelStim.stmhTqHist); 
rStartRelStim.stmhTqHcNormCumSum = cumsum(rStartRelStim.stmhTqHistNorm); % cumulative sum
bar(edges(1:end-1),rStartRelStim.stmhTqHistNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
set(gca,'TickDir','out')

% high torque pStim
[rStartRelStim.pStmhTqHist, rStartRelStim.pStmhTqHistE] = histcounts(rStartRelStim.pStmhTq, edges);
rStartRelStim.pStmhTqHistNorm = rStartRelStim.pStmhTqHist./sum(rStartRelStim.pStmhTqHist); 
rStartRelStim.pStmhTqHcNormCumSum = cumsum(rStartRelStim.pStmhTqHistNorm); % cumulative sum
bar(edges(1:end-1),rStartRelStim.pStmhTqHistNorm,1,'FaceAlpha',.5,'EdgeColor','none'); 
hold off; 
ylim([0 max(rStartRelStim.pStmhTqHistNorm)+.02])
yticks(0:.1:1);
set(gca,'TickDir','out')

% plot cumulative pullStart on high torque
figure; hold on; 
plot(edges(1:end-1), rStartRelStim.pStmhTqHcNormCumSum)
plot(edges(1:end-1), rStartRelStim.stmhTqHcNormCumSum)
set(gca,'TickDir','out')

end






















