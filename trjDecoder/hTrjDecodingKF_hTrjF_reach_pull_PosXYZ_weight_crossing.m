function hTrjDecodingKF_hTrjF_reach_pull_PosXYZ_weight_crossing(filePath)
%This decodes kinematics of mouse 3-d hand movement trajectories (X,Y,Z)
% using cross-validated (leave-a-trial-out) Kalman filter decoding.

cd(filePath)
kfDir = dir('preprocessKFdecodeHTrj_reachpull_new*');

if ~isempty(kfDir)
    
    wrI = strfind(filePath, 'WR');
    m_name = filePath(wrI:wrI+10);
    
    load(fullfile(kfDir.folder,kfDir.name),'s')
    
    ctx_record = isfield(s.dat, 'spkCtxR');
    str_record = isfield(s.dat, 'spkStrR');
    cg_record = isfield(s.dat, 'spkCgR');
    
    resample = 1;
    
    if ctx_record && str_record
        valTrI_reach = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxR, 'un', 0));
        valTrI_pull = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxP, 'un', 0));
        valTrI = valTrI_reach & valTrI_pull;
    elseif cg_record  % cg only
        valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCgR, 'un', 0));
    end
    
    %% REACH select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
    tmpState0R = cell(size(s.dat.stateR,1),size(s.dat.stateR,2));
    tmpStateR = cell(size(s.dat.stateR,1),size(s.dat.stateR,2));
    
    tmpState0P = cell(size(s.dat.stateP,1),size(s.dat.stateP,2));
    tmpStateP = cell(size(s.dat.stateP,1),size(s.dat.stateP,2));
    
    tmpState0R(valTrI)=deal(cellfun(@(a) a(1:3,:),s.dat.stateR(valTrI), 'un', 0));
    tmpState0P(valTrI)=deal(cellfun(@(a) a(1:3,:),s.dat.stateP(valTrI), 'un', 0));
    
    % global median subtraction (not baseline subtraction of its own)
    medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), tmpState0R(valTrI)', 'un', 0)),2);
    tmpStateR(valTrI) = cellfun(@(a) a-repmat(medP1, 1, size(a,2)), tmpState0R(valTrI), 'un', 0);
    tmpStateP(valTrI) = cellfun(@(a) a-repmat(medP1, 1, size(a,2)), tmpState0P(valTrI), 'un', 0);

    %% leave-a-trial-out decoding using Kalman Filter (heavy-lifting part)
    rez = leaveOneOutKFdecoderTrialTypeBalanced_withCg_weight_crossing(tmpStateR, tmpStateP, s);
    
    %% evaluate decoding with correlation and r-squred
    if ctx_record && str_record && cg_record  % all three regions
        [corrRez.ctx, r2Rez.ctx] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateCtxM);
        [corrRez.str, r2Rez.str] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateStrM);
        [corrRez.cg, r2Rez.cg] = trjDecodeEvalByTrialType(tmpStateP, rez.estStateCgM);
        
        [corrRez.ctx_mismatch, r2Rez.ctx_mismatch] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateCtxM_mismatch);
        [corrRez.str_mismatch, r2Rez.str_mismatch] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateStrM_mismatch);
        [corrRez.cg_mismatch, r2Rez.cg_mismatch] = trjDecodeEvalByTrialType(tmpStateP, rez.estStateCgM_mismatch);
    elseif ctx_record && str_record && ~cg_record  % ctx & str
        [corrRez.ctx, r2Rez.ctx] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateCtxM);
        [corrRez.str, r2Rez.str] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateStrM);
        
        [corrRez.ctx_mismatch, r2Rez.ctx_mismatch] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateCtxM_mismatch);
        [corrRez.str_mismatch, r2Rez.str_mismatch] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateStrM_mismatch);
    elseif ~ctx_record && ~str_record && cg_record  % cg only
        [corrRez.cg, r2Rez.cg] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateCgM);
        [corrRez.cg_mismatch, r2Rez.cg_mismatch] =  trjDecodeEvalByTrialType(tmpStateP, rez.estStateCgM_mismatch);
    end
    
%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_weight_mismatch_', m_name)),'s', 'rez', 'corrRez','r2Rez')     
    
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix interpolation function
function [intMat] = intm(origMat, numbDataPoints)
% 1-d interpolation of a matrix as specified by the number of data points
% numbDataPoints = 100; % 20ms*100 = 2000ms
x=1:size(origMat,2);
xq=linspace(1,size(origMat,2),numbDataPoints);
intMat = interp1(x,origMat',xq)';
if size(intMat,1)>size(intMat,2)
    intMat = intMat'; % to ensure variable-by-time orientation
end
end

function [corrRez,r2Rez] = trjDecodeEvalByTrialType(state, stateEst)
%Computes pearson correlation and r-squared metrics between actual and
% estimated states, each input should be trial-by-trialType cell, whose cell
% element should be variables-by-timeBin.
% state = s.dat.stateP; % actual state
% stateEst = s.dat.statePCtx; % estimated state

nEach=sum(cell2mat(cellfun(@(a) ~isempty(a), state, 'un',0)),1);
valType = nEach>10; % evalute trial types with more than 10 trials

isNaN_stateEst = cell2mat(cellfun(@(a) sum(sum(isnan(a)))>=1, stateEst, 'un', 0));
if sum(isNaN_stateEst(:)) < 10
    
    if sum(isNaN_stateEst(:)) > 0  % to deal with NaN values
        [stateEst{isNaN_stateEst}] = deal([]);
        [state{isNaN_stateEst}] = deal([]);
    end
    
    state1 = state(:,valType);
    stateEst1 = stateEst(:,valType);
    
    s1 = cell2mat(reshape(state1,[],1)')'; % state mat
    s2 = cell2mat(reshape(stateEst1,[],1)')'; % stateEst mat
    
    r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
    r2overall = @(a,b) 1-sum(nansum((a-b).^2))./sum(nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2)); % r-squared = 1-SSres/SStot;
    r2Rez.all = r2(s1,s2); % r^2 across all trials
    r2Rez.overall = r2overall(s1,s2);
    
    corrRez.all = corr(s1,s2,'Rows','complete'); % corr without smoothing
    corrRez.all_sm = corr(smooth2a(s1,3,0),smooth2a(s2,3,0),'Rows','complete'); % corr with smoothing across time bins (smoothing seems to hurt usually)
    
    for tt = 1:size(state,2) % trial types
        if tt<= size(state1,2)
            s1tt = cell2mat(reshape(state(:,tt),[],1)')';
            s2tt = cell2mat(reshape(stateEst(:,tt),[],1)')';
            tempCorr = corr(s1tt,s2tt,'Rows','complete');
            tempCorr_sm = corr(smooth2a(s1tt,3,0),smooth2a(s2tt,3,0),'Rows','complete');
            tempR2 = r2(s1tt,s2tt);
        else
            tempCorr = nan;
            tempCorr_sm = nan;
            tempR2 = nan;
        end
        
        switch tt
            case 1
                corrRez.lelt = tempCorr;
                corrRez.lelt_sm = tempCorr_sm;
                r2Rez.lelt = tempR2;
            case 2
                corrRez.leht = tempCorr;
                corrRez.leht_sm = tempCorr_sm;
                r2Rez.leht = tempR2;
            case 3
                corrRez.rilt = tempCorr;
                corrRez.rilt_sm = tempCorr_sm;
                r2Rez.rilt = tempR2;
            case 4
                corrRez.riht = tempCorr;
                corrRez.riht_sm = tempCorr_sm;
                r2Rez.riht = tempR2;
        end
    end
else  % in case there are too many NaNs
    corrRez = NaN;
    r2Rez = NaN;
    
end
end


















