function reevaluateDecodeRezByTrialType(filePath, saveName)
%This decodes kinematics of mouse 3-d hand movement trajectories (X,Y,Z)
% using cross-validated (leave-a-trial-out) Kalman filter decoding.
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)
kfDir = dir('rezKFdecodeHTrjCtxStrPosVel_reach_WR*');
load(fullfile(kfDir.folder,kfDir.name),'s')
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxR, 'un', 0));
tmpState0 = cell(size(s.dat.stateR,1),size(s.dat.stateR,2)); 
tmpState = cell(size(s.dat.stateR,1),size(s.dat.stateR,2)); 
for k = 1:8 % MAIN LOOP (fit X,Y,Z sperately and XYZ altogether) 
    if k <= 6
        tmpState0(valTrI)=deal(cellfun(@(a) a(k,:),s.dat.stateR(valTrI),'un',0)); 
    elseif k == 7 % position variables together
        tmpState0(valTrI)=deal(cellfun(@(a) a(1:3,:),s.dat.stateR(valTrI),'un',0)); 
    elseif k == 8 % velosity variables together
        tmpState0(valTrI)=deal(cellfun(@(a) a(4:6,:),s.dat.stateR(valTrI),'un',0)); 
    end
    medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), tmpState0(valTrI)','un',0)),2);
    tmpState(valTrI) = cellfun(@(a) a-repmat(medP1,1,size(a,2)),tmpState0(valTrI),'un',0);      
    %% evaluate decoding with correlation and r-squred 
    [corrRez.ctx{k},r2Rez.ctx{k}] =  trjDecodeEvalByTrialType( tmpState, s.dat.stateRCtx{k}.est); 
    [corrRez.str{k},r2Rez.str{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.stateRStr{k}.est); 
    [corrRez.ctxstr{k},r2Rez.ctxstr{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.stateRCtxStr{k}.est); 
end
clearvars k

%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_reach_',saveName)),'corrRez','r2Rez') % last saved after training without stim trials 5/27 Wed 9pm
%trjMovie([stateCtxCC_sm(:,2), stateStrCC_sm(:,2), stateCC_sm(:,2)]', figSaveDir, 'kfDecode_Ypos_CtxStrAct')

%% PULL
kfDir = dir('rezKFdecodeHTrjCtxStrPosVel_pull_WR*');
load(fullfile(kfDir.folder,kfDir.name),'s')
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxP, 'un', 0));

%% REACH select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
tmpState0 = cell(size(s.dat.stateP,1),size(s.dat.stateP,2)); 
tmpState = cell(size(s.dat.stateP,1),size(s.dat.stateP,2)); 
for k = 1:8 % MAIN LOOP (fit X,Y,Z sperately and XYZ altogether) 
    if k <= 6
        tmpState0(valTrI)=deal(cellfun(@(a) a(k,:),s.dat.stateP(valTrI),'un',0)); 
    elseif k == 7 % position variables together
        tmpState0(valTrI)=deal(cellfun(@(a) a(1:3,:),s.dat.stateP(valTrI),'un',0)); 
    elseif k == 8 % velosity variables together
        tmpState0(valTrI)=deal(cellfun(@(a) a(4:6,:),s.dat.stateP(valTrI),'un',0)); 
    end
    medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), tmpState0(valTrI)','un',0)),2);
    tmpState(valTrI) = cellfun(@(a) a-repmat(medP1,1,size(a,2)),tmpState0(valTrI),'un',0);      
    %% evaluate decoding with correlation and r-squred 
    [corrRez.ctx{k},r2Rez.ctx{k}] =  trjDecodeEvalByTrialType( tmpState, s.dat.statePCtx{k}.est); 
    [corrRez.str{k},r2Rez.str{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.statePStr{k}.est); 
    [corrRez.ctxstr{k},r2Rez.ctxstr{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.statePCtxStr{k}.est); 
end
clearvars k

%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_pull_',saveName)),'corrRez','r2Rez') % last saved after training without stim trials 5/27 Wed 9pm
