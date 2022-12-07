function m1StimEffectAnalysisJsTrj(filePath)

if contains(filePath, 'BehVariablesJs')
    S = load(fullfile(filePath), 'jsTime1k');  % load jsTime1k_KV.mat
    S = S.('jsTime1k');
elseif contains(filePath, 'jsTime1k_Kinematics_VideoFiles')
    S = load(fullfile(filePath), 'jsTime1k_KV');  % load jsTime1k_KV.mat
    S = S.('jsTime1k_KV');
elseif contains(filePath, 'jsTime1k_KinematicsTrajectories')
    S = load(fullfile(filePath), 'jkvt');  % load jsTime1k_KV.mat
    S = S.('jkvt');
end

firstRwdTrial = find([S(:).rewarded]==1,1,'first'); 

m_name = filePath(strfind(filePath, 'WR'):strfind(filePath, 'WR')+10); 

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
            plot(1:trStartPt+length(S(t).movKins.sgJsTrajmm)-trStartPt, S(t).movKins.sgJsTrajmm, 'Color',[70,240,240]./255);
            plot(trStartPt+1-trStartPt, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [70,240,240]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
            %plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[70,240,240]./255);
            %plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [70,240,240]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end
    elseif ~isnan(S(t).stimLaserOn) % a stim trial
        if isstruct(S(t).movKins)
            trStartPt = S(t).trJsReady-S(t).stimLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            plot(1:length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[240,50,230]./255);
            plot(1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [240,50,230]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);            
            %plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[240,50,230]./255);
            %plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [240,50,230]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end     
    end
    
end
clearvars t trStartPt
hold off;
xlim([-1000 5000]); xlabel('Time relative to actual/putative stim onset')
ylim([-20 20]); ylabel('Joystick position (mm)')
set(gca, 'TickDir', 'out')
print(fullfile(fileparts(filePath),'Figure','JsTrajStimVsNoStim_AllTrials'), '-dpdf','-painters')

% Successful trials only
hold on; 
for t = 2:length(S)
    
    if isnan(S(t).stimLaserOn)
        if isstruct(S(t).movKins) && strcmpi(S(t).trialType,'sp')
            if isfield(S,'pLaserOn') % if pseudo laser pulse was used 
              %trStartPt = S(t).trJsReady-S(t).pLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
              %trStartPtCollect(t,1) = trStartPt; 
            else
              %trStartPt = S(t).trJsReady-(S(t-1).trEnd+4000); % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            end
            plot(1:length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[70,240,240]./255);
            plot(1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [70,240,240]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);   
            %plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[70,240,240]./255);
            %plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [70,240,240]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end
    elseif ~isnan(S(t).stimLaserOn) && strcmpi(S(t).trialType,'sp')
        if isstruct(S(t).movKins)
            %trStartPt = S(t).trJsReady-S(t).stimLaserOn; % align to the putative stim on, stim generated 4-sec after the last trial offset with the camTrigger
            plot(1:length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[240,50,230]./255);
            plot(1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [240,50,230]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
            %plot(trStartPt+1:trStartPt+length(S(t).movKins.sgJsTrajmm), S(t).movKins.sgJsTrajmm, 'Color',[240,50,230]./255);
            %plot(trStartPt+1, S(t).movKins.sgJsTrajmm(1), 'o', 'MarkerFaceColor', [240,50,230]./255, 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end     
    end
   
end
clearvars t trStartPt
xlim([-1000 5000]); xlabel('Time relative to actual/putative stim onset(ms)')
ylim([-20 5]); ylabel('Joystick position (mm)')
set(gca, 'TickDir', 'out')
print(fullfile(fileparts(filePath),'Figure','JsTrajStimVsNoStim-SuccessTrials'), '-dpdf','-painters')
hold off;

%% Reach probability P(reachOn|Time), Expected Js position E(JsPosition|Time)
valTrCnt = 0; 
for t = 1:length(S)
  
    if isstruct(S(t).movKins) && t>firstRwdTrial
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

g(1,2).stat_summary();
g(1,2).set_title('P(reachOn|Time)');

figure('Position',[100 100 800 550]);
g.draw();
%set(gcf,'renderer','painters')
print(fullfile(fileparts(filePath),'Figure','reachProb_expectedJsPos'), '-dpdf','-bestfit','-painters')

stimI = cell2mat(cellfun(@(a) strcmpi(a, 'stim'), stimTrIdx, 'un', 0)); 
jsReachOn_stim = cell2mat(jsReachOn(stimI)');
jsReachOn_prob_stim = sum(sum(jsReachOn_stim(:, 1:1000), 2) > 0)/size(jsReachOn_stim, 1); 

jsReachOn_noStim = cell2mat(jsReachOn(~stimI)'); 
jsReachOn_prob_noStim = sum(sum(jsReachOn_noStim(:, 1:1000), 2) > 0)/size(jsReachOn_noStim, 1); 

save(fullfile(fileparts(filePath), strcat('stim_effect_rez_', m_name)), 'jsPosE', 'stimTrIdx', 'jsReachOn', 'jsReachOn_prob_stim', 'jsReachOn_prob_noStim')

end


