function hTrjDecodingKalmanFilter_hTrjF_reach_pull_PosVelXYZ_params(filePath, saveName)
%This decodes kinematics of mouse 3-d hand movement trajectories (X,Y,Z)
% using cross-validated (leave-a-trial-out) Kalman filter decoding.
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)
kfDir = dir('preprocessKFdecodeHTrjCtxStr_reachpull_hTrjF_20ms*');
load(fullfile(kfDir.folder,kfDir.name),'s')

resample = 50; 
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxR, 'un', 0));

%% REACH select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
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
    % global median subtraction (not baseline subtraction of its own)
    medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), tmpState0(valTrI)','un',0)),2);
    tmpState(valTrI) = cellfun(@(a) a-repmat(medP1,1,size(a,2)),tmpState0(valTrI),'un',0);
    %% leave-a-trial-out decoding using Kalman Filter (heavy-lifting part)
    [s.dat.stateRCtx{k},s.dat.stateRStr{k},s.dat.stateRCtxStr{k},s.dat.params{k}] = leaveOneOutKFdecoderTrialTypeBalanced_params(tmpState, s.dat.spkCtxR, s.dat.spkStrR, s.dat.laserIdxR, resample, s.ctxI, s.strI);
    
    %% evaluate decoding with correlation and r-squred
    [corrRez.ctx{k},r2Rez.ctx{k}] =  trjDecodeEvalByTrialType(tmpState, s.dat.stateRCtx{k}.est);
    [corrRez.str{k},r2Rez.str{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.stateRStr{k}.est);
    [corrRez.ctxstr{k},r2Rez.ctxstr{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.stateRCtxStr{k}.est);
end
clearvars k

%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_reach_',saveName)),'s','corrRez','r2Rez') % last saved after training without stim trials 5/27 Wed 9pm
%trjMovie([stateCtxCC_sm(:,2), stateStrCC_sm(:,2), stateCC_sm(:,2)]', figSaveDir, 'kfDecode_Ypos_CtxStrAct')

%% PULL
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxP, 'un', 0));

%% PULL select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
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
    % global median subtraction (not baseline subtraction of its own)
    medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), tmpState0(valTrI)','un',0)),2);
    tmpState(valTrI) = cellfun(@(a) a-repmat(medP1,1,size(a,2)),tmpState0(valTrI),'un',0);
    %% leave-a-trial-out decoding using Kalman Filter (heavy-lifting part)
    %[s.dat.stateRCtx{k},s.dat.stateRStr{k},s.dat.stateRCtxStr{k}] = leaveOneOutKFdecoder(tmpState, s.dat.spkCtxP, s.dat.spkStrP, s.dat.laserIdxP, 10);
    [s.dat.stateRCtx{k},s.dat.stateRStr{k},s.dat.stateRCtxStr{k},s.dat.params{k}] = leaveOneOutKFdecoderTrialTypeBalanced_params(tmpState, s.dat.spkCtxP, s.dat.spkStrP, s.dat.laserIdxP, resample, s.ctxI, s.strI);
    
    %% evaluate decoding with correlation and r-squred
    [corrRez.ctx{k},r2Rez.ctx{k}] =  trjDecodeEvalByTrialType(tmpState, s.dat.stateRCtx{k}.est);
    [corrRez.str{k},r2Rez.str{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.stateRStr{k}.est);
    [corrRez.ctxstr{k},r2Rez.ctxstr{k}] = trjDecodeEvalByTrialType(tmpState, s.dat.stateRCtxStr{k}.est);
end
clearvars k

%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_pull_',saveName)),'s','corrRez','r2Rez') % last saved after training without stim trials 5/27 Wed 9pm

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
    end

%     function [corrRez,r2Rez] = trjDecodeEvalByTrialType(state, stateEst)
%         %Computes pearson correlation and r-square metrics between actual and
%         % estimated states, each input should be trial-by-trialType cell, whose cell
%         % element should be variables-by-timeBin.
%         % state = s.dat.stateP; % actual state
%         % stateEst = s.dat.statePCtx; % estimated state
%         s1 = cell2mat(reshape(state,[],1)')'; % state mat
%         s2 = cell2mat(reshape(stateEst,[],1)')'; % stateEst mat
%         
%         r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
%         r2Rez.all = r2(s1,s2); % r^2 across all trials
%         
%         corrRez.all = corr(s1,s2,'Rows','complete'); % corr without smoothing
%         corrRez.all_sm = corr(smooth2a(s1,3,0),smooth2a(s2,3,0),'Rows','complete'); % corr with smoothing across time bins (smoothing seems to hurt usually)
%         
%         for tt = 1:size(state,2) % trial types
%             s1tt = cell2mat(reshape(state(:,tt),[],1)')';
%             s2tt = cell2mat(reshape(stateEst(:,tt),[],1)')';
%             tempCorr = corr(s1tt,s2tt,'Rows','complete');
%             tempCorr_sm = corr(smooth2a(s1tt,3,0),smooth2a(s2tt,3,0),'Rows','complete');
%             tempR2 = r2(s1tt,s2tt);
%             
%             switch tt
%                 case 1
%                     corrRez.lelt = tempCorr;
%                     corrRez.lelt_sm = tempCorr_sm;
%                     r2Rez.lelt = tempR2;
%                 case 2
%                     corrRez.leht = tempCorr;
%                     corrRez.leht_sm = tempCorr_sm;
%                     r2Rez.leht = tempR2;
%                 case 3
%                     corrRez.rilt = tempCorr;
%                     corrRez.rilt_sm = tempCorr_sm;
%                     r2Rez.rilt = tempR2;
%                 case 4
%                     corrRez.riht = tempCorr;
%                     corrRez.riht_sm = tempCorr_sm;
%                     r2Rez.riht = tempR2;
%             end
%         end
%     end
end















