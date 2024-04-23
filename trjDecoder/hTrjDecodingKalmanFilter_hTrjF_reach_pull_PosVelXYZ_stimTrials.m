function hTrjDecodingKalmanFilter_hTrjF_reach_pull_PosVelXYZ_stimTrials(filePath, saveName)
%This decodes kinematics of mouse 3-d hand movement trajectories (X,Y,Z)
% using cross-validated (leave-a-trial-out) Kalman filter decoding.
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles';
cd(filePath)

kfDir = dir('preprocessKFdecodeHTrj_reachpull_new*');

if ~isempty(kfDir)
    % load decoding result for reach phase (rez_reach)
    load(fullfile(filePath, strcat('rezKFdecodeHTrjCtxStrPosVel_reach_new_', saveName)), 's', 'rez_reach')

    % load decoding result for pull phase (rez_pull)
    load(fullfile(filePath, strcat('rezKFdecodeHTrjCtxStrPosVel_pull_new_', saveName)), 'rez_pull')

    ctx_logic = isfield(s.dat, 'spkCtxR');
    str_logic = isfield(s.dat, 'spkStrR');
    cg_logic = isfield(s.dat, 'spkCgR');

    %% get valTrI
    if ctx_logic && str_logic
        valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxR, 'un', 0));
    elseif cg_logic  % cg only
        valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCgR, 'un', 0));
    end

    %% get state (position) during reach
    % raw reach trace
    statePosRawRch = cell(size(s.dat.stateR,1),size(s.dat.stateR,2));
    statePosRawRch(valTrI)=deal(cellfun(@(a) a(1:3,:),s.dat.stateR(valTrI),'un',0));

    % median initial reach position subtracted reach trace
    statePosRch = cell(size(s.dat.stateR,1),size(s.dat.stateR,2));
    medPos1Rch = nanmedian(cell2mat(cellfun(@(a) a(:,1), statePosRawRch(valTrI)','un',0)), 2);
    statePosRch(valTrI) = cellfun(@(a) a-repmat(medPos1Rch, 1, size(a, 2)), statePosRawRch(valTrI),'un',0);

    % raw pull trace
    statePosRawPull = cell(size(s.dat.stateP,1),size(s.dat.stateP,2));
    statePosRawPull(valTrI)=deal(cellfun(@(a) a(1:3,:),s.dat.stateP(valTrI),'un',0));

    % median initial pull position subtracted pull trace
    statePosPull = cell(size(s.dat.stateP,1),size(s.dat.stateP,2));
    medPos1Pull = nanmedian(cell2mat(cellfun(@(a) a(:,1), statePosRawPull(valTrI)','un',0)), 2);
    statePosPull(valTrI) = cellfun(@(a) a-repmat(medPos1Pull, 1, size(a, 2)), statePosRawPull(valTrI),'un',0);

    evaluateStimTrialDecoding(statePosRch, rez_reach{7}.estStateCtxM, s.dat.laserIdxR)

    evaluateStimTrialDecoding(statePosPull, rez_pull{7}.estStateCtxM, s.dat.laserIdxP)

    function evaluateStimTrialDecoding(stateRch, statePull, estStateRch, estStatePull, mPos1Rch, mPos1Pull, stimRchI, stimPullI)
%             stateRch = statePosRch;
%             statePull = statePosPull;
%             estStateRch = rez_reach{7}.estStateCtxM;
%             estStatePull = rez_pull{7}.estStateCtxM;
%             stimRchI = s.dat.laserIdxR;
%             stimPullI = s.dat.laserIdxP;
%             mPos1Rch = medPos1Rch;
%             mPos1Pull = medPos1Pull;
        %     resampled distribution
        %     Q1. Does silencing slows the speed of reach and/or pull trajectory?

        %     define the starting point 'p0' based on the prior knowledge (100 ms baseline period, 5 timebins), position
        %     and velocity

        %     get expected errors for each trajectory 'steps'
        %       1) Define the starting point and reach end point
        %       2) Use interpolation to get equal numbers of steps per trajectory


        % stim trial logic
        stimTrI = cellfun(@sum, stimRchI)>0 | cellfun(@sum, stimPullI)>0;

        % XYZ alignment to be added to the pull trajectories to recombine the reach and pull trajectories
        alignRchPull = mPos1Pull-mPos1Rch;

        % interpolated stim logic
        stimRchIntI = cell(size(stimRchI));
        stimPullIntI = cell(size(stimPullI));

        statePullAlign = cell(size(statePull)); % pull trj aligned to median initial reach point to combine reach and pull trajectories
        stateRchPull = cell(size(stateRch)); % for combined reach and pull trajectories
        stateRchInt = cell(size(stateRch)); % interpolated reach trajectories
        statePullInt = cell(size(statePull)); % interpolated pull trajectories

        statePullAlignEst = cell(size(estStatePull)); % estimated pull trj aligned to median intial reach point
        stateRchPullEst = cell(size(estStateRch)); % for combined reach and pull estimated trajectories
        stateRchIntEst = cell(size(estStateRch)); % interpolated estimated reach trajectories
        statePullIntEst = cell(size(estStatePull)); % interpolated estimated pull trajectories

        distRchEst = cell(size(estStateRch)); % Euclidean distance of estimated reach trajectories
        distPullEst = cell(size(estStatePull)); % Euclidean distance of estimated pull trajectories

        for r = 1:size(stateRch, 1)
            for c = 1:size(stateRch, 2)
                if ~isempty(stateRch{r, c})
                    rTrj = stateRch{r, c};
                    rTrjEst = estStateRch(r, c, :); % resampled

                    pTrj = statePull{r, c};
                    pTrjEst = estStatePull(r, c, :); % resampled

                    % recombine the reach and pull trajectories
                    statePullAlign{r, c} = pTrj + repmat(alignRchPull, 1, size(pTrj, 2));

                    % interpolate reach start to end
                    [~, r0I] = min(sum(rTrj.^2));
                    if length(rTrj)-r0I>=3
                        stateRchInt{r, c} = intm(rTrj(:, r0I:end), 10);
                        stateRchIntEst(r, c, :) = cellfun(@(a) intm(a(:, r0I:end), 10), rTrjEst, 'UniformOutput', false);
                        distRchEst(r, c, :) = cellfun(@(a) sqrt((a-stateRchInt{r, c}).^2), stateRchIntEst(r, c, :), 'UniformOutput', false);

                        if stimTrI(r, c)
                            stimRchIntI{r, c} = intm(double(stimRchI{r, c}(r0I:end)), 10)>0;
                        end
                    end

                    % interpolate pull start to end
                    [~, p0I] = min(sum(pTrj.^2));
                    if length(pTrj)-p0I>=3
                        statePullInt{r, c} = intm(pTrj(:, p0I:end), 10);
                        statePullIntEst(r, c, :) = cellfun(@(a) intm(a(:, p0I:end), 10), pTrjEst, 'UniformOutput', false);
                        distPullEst(r, c, :) = cellfun(@(a) sqrt((a-statePullInt{r, c}).^2), statePullIntEst(r, c, :), 'UniformOutput', false);

                        if stimTrI(r, c)
                            stimPullIntI{r, c} = intm(double(stimPullI{r, c}(p0I:end)), 10)>0;
                        end
                    end
                end
            end
        end

        %%
        stimTrI3 = repmat(stimTrI, 1, 1, size(distRchEst, 3));

        % no-stim reach
        distIntRchEstNoStimC = distRchEst(~stimTrI3);
        distIntRchEstNoStim = [];

        for i = 1:length(distIntRchEstNoStimC)
            if ~isempty(distIntRchEstNoStimC{i})
                distIntRchEstNoStim(:, :, end+1) = distIntRchEstNoStimC{i};
            end
        end
        mDistIntRchEstNoStim = mean(distIntRchEstNoStim, 3);

        % stim reach
        distIntRchEstStimC = distRchEst(stimTrI3);
        distIntRchEstStim = [];

        for ii = 1:length(distIntRchEstStimC)
            if ~isempty(distIntRchEstStimC{ii})
                distIntRchEstStim(:, :, end+1) = distIntRchEstStimC{ii};
            end
        end
        mDistIntRchEstStim = mean(distIntRchEstStim, 3);

        %mRelDistRchEst = cell(size(stimTrI)); 
        mRelDistRchEst = nan(3, 10, sum(stimTrI(:))); 
        countStimTr = 0; 
        % stim reach points
        for r = 1:size(distRchEst, 1)
            for c = 1:size(distRchEst, 2)
                if stimTrI(r, c) && sum(~cellfun(@isempty, distRchEst(r, c, :)))==length(distRchEst(r, c, :))
                    countStimTr = countStimTr + 1; 
                    mRelDistRchEst(:, :, countStimTr) = mean(cell2mat(cellfun(@(a) a./mDistIntRchEstNoStim, distRchEst(r, c, :), 'UniformOutput', false)), 3);
                    mRelDistRchEst(:, ~stimRchIntI{r, c}, countStimTr) = NaN;
                end
            end
        end

        nanmean(mRelDistRchEst, 3)

        % no-stim pull
        distIntPullEstNoStimC = distRchEst(~stimTrI3);
        distIntPullEstNoStim = [];

        for j = 1:length(distIntPullEstNoStimC)
            if ~isempty(distIntPullEstNoStimC{j})
                distIntPullEstNoStim(:, :, end+1) = distIntPullEstNoStimC{j};
            end
        end
        mDistIntPullEstNoStim = mean(distIntPullEstNoStim, 3);

        % stim pull
        distIntPullEstStimC = distRchEst(stimTrI3);
        distIntPullEstStim = [];

        for jj = 1:length(distIntPullEstStimC)
            if ~isempty(distIntPullEstStimC{jj})
                distIntPullEstStim(:, :, end+1) = distIntPullEstStimC{jj};
            end
        end
        mDistIntPullEstStim = mean(distIntPullEstNoStim, 3);

        mRelDistPullEst = nan(3, 10, sum(stimTrI(:))); 
        countStimTr = 0; 
        % stim pull points
        for r = 1:size(distPullEst, 1)
            for c = 1:size(distPullEst, 2)
                if stimTrI(r, c) && sum(~cellfun(@isempty, distPullEst(r, c, :)))==length(distPullEst(r, c, :))
                    countStimTr = countStimTr + 1; 
                    mRelDistPullEst(:, :, countStimTr) = mean(cell2mat(cellfun(@(a) a./mDistIntPullEstNoStim, distPullEst(r, c, :), 'UniformOutput', false)), 3);
                    mRelDistPullEst(:, ~stimPullIntI{r, c}, countStimTr) = NaN;
                end
            end
        end

        nanmean(mRelDistPullEst, 3)













    end
end


    function error = trjDistErr(trj, trjEst)
        %trj = rTrj;
        %trjEst = rTrjEst;

        % Calculate the squared differences
        sqDiffs = (trj - trjEst).^2;

        % Take the square root to get the Euclidean distance
        euD = sqrt(sqDiffs);

        % Compute basic statistics
        error.mean = mean(euD, 2);
        error.median = median(euD, 2);
        error.std = std(euD, 0, 2);
    end


end





%% evaluate decoding with correlation and r-squared
if ctx_record && str_record && cg_record  % all three regions
    [corrRez_reach.ctx{k}, r2Rez_reach.ctx{k}] =  trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_ctx);
    [corrRez_reach.ctx1{k}, r2Rez_reach.ctx1{k}] =  trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_ctx1);  % cell # matched to Cg
    [corrRez_reach.str{k}, r2Rez_reach.str{k}] =  trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_str);
    [corrRez_reach.str1{k}, r2Rez_reach.str1{k}] =  trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_str1);  % cell # matched to Cg
    [corrRez_reach.cg{k}, r2Rez_reach.cg{k}] = trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_cg);
    [corrRez_reach.ctx_str{k}, r2Rez_reach.ctx_str{k}] = trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_ctx_str);
    [corrRez_reach.ctx_str_cg{k}, r2Rez_reach.ctx_str_cg{k}] = trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_ctx_str_cg);
elseif ctx_record && str_record && ~cg_record  % ctx & str
    [corrRez_reach.ctx{k}, r2Rez_reach.ctx{k}] =  trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_ctx);
    [corrRez_reach.str{k}, r2Rez_reach.str{k}] =  trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_str);
    [corrRez_reach.ctx_str{k}, r2Rez_reach.ctx_str{k}] = trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_ctx_str);
elseif ~ctx_record && ~str_record && cg_record  % cg only
    [corrRez_reach.cg{k}, r2Rez_reach.cg{k}] = trjDecodeEvalByTrialType(tmpState, rez_reach{k}.rpr_est_cg);
end
end
clearvars k

%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_reach_new_', saveName)), 's', 'rez_reach', 'corrRez_reach', 'r2Rez_reach') % last saved after training without stim trials 5/27 Wed 9pm
%trjMovie([stateCtxCC_sm(:,2), stateStrCC_sm(:,2), stateCC_sm(:,2)]', figSaveDir, 'kfDecode_Ypos_CtxStrAct')

%% PULL
if ctx_record && str_record
    valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtxP, 'un', 0));
elseif cg_record  % cg only
    valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCgP, 'un', 0));
end

%% PULL select kinematic variables to fit (e.g. hand position or hand velocity - fitting them both together doesn't seem to be a good idea for some reason(?))
tmpState0 = cell(size(s.dat.stateP,1),size(s.dat.stateP,2));
tmpState = cell(size(s.dat.stateP,1),size(s.dat.stateP,2));
for k = 7:8 % MAIN LOOP (fit X,Y,Z sperately and XYZ altogether)
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
    rez_pull{k} = leaveOneOutKFdecoderTrialTypeBalanced_pull_withCg(tmpState, s, resample);

    %% evaluate decoding with correlation and r-squred
    if ctx_record && str_record && cg_record  % all three regions
        [corrRez_pull.ctx{k}, r2Rez_pull.ctx{k}] =  trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_ctx);
        [corrRez_pull.ctx1{k}, r2Rez_pull.ctx1{k}] =  trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_ctx1);  % cell # matched to Cg
        [corrRez_pull.str{k}, r2Rez_pull.str{k}] =  trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_str);
        [corrRez_pull.str1{k}, r2Rez_pull.str1{k}] =  trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_str1);  % cell # matched to Cg
        [corrRez_pull.cg{k}, r2Rez_pull.cg{k}] = trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_cg);
        [corrRez_pull.ctx_str{k}, r2Rez_pull.ctx_str{k}] = trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_ctx_str);
        [corrRez_pull.ctx_str_cg{k}, r2Rez_pull.ctx_str_cg{k}] = trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_ctx_str_cg);
    elseif ctx_record && str_record && ~cg_record  % ctx & str
        [corrRez_pull.ctx{k}, r2Rez_pull.ctx{k}] =  trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_ctx);
        [corrRez_pull.str{k}, r2Rez_pull.str{k}] =  trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_str);
        [corrRez_pull.ctx_str{k}, r2Rez_pull.ctx_str{k}] = trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_ctx_str);
    elseif ~ctx_record && ~str_record && cg_record  % cg only
        [corrRez_pull.cg{k}, r2Rez_pull.cg{k}] = trjDecodeEvalByTrialType(tmpState, rez_pull{k}.rpr_est_cg);
    end
end
clearvars k

%% save the result
save(fullfile(filePath,strcat('rezKFdecodeHTrjCtxStrPosVel_pull_new_',m_name)),'s', 'rez_pull', 'corrRez_pull','r2Rez_pull') % last saved after training without stim trials 5/27 Wed 9pm

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


end
















