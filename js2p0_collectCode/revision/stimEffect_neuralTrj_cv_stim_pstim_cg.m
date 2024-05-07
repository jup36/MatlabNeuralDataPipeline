function [rrrCv, rrrRezR2] = stimEffect_neuralTrj_cv_stim_pstim_cg(filePath, dims, folds)
%This function is to run the reduced rank regression (RRR) with cross-validation,
% which trains the RRR model that predicts the target (e.g., str) neural population activity
% based on the source neural population activity (e.g., ctx) with cross validation.
% The trained models are tested both on the held-out no-stim and stim
% trials. The results of rrr and testing (r2) are saved and/or returned as
% an output of this function.

saveNameIdx = strfind(filePath, 'WR');
saveName = filePath(saveNameIdx(1):saveNameIdx(1)+10);
saveDir = fileparts(filePath);

%filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/js2p0_tbytSpkHandJsTrjBin_50ms_stimPstimPrepExtWoTo_WR40_081919.mat';
load(fullfile(filePath), 'ss')

%% get indices
% valid reach and pull trials (trials with valid hand trajectories)
valRchTrI = cellfun(@(a) ~isempty(a), {ss.spkRchIdx});
valPullTrI = cellfun(@(a) ~isempty(a), {ss.spkPullIdx});

excludeNaNsInBothPops = @(a, b) sum(sum(isnan(full(a)), 2))==0 & sum(sum(isnan(full(b)), 2))==0;

% stim trials
stimTrI =  cellfun(@(a) ~isempty(a), {ss.utbCgStimAlign}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCgStimAlign}, {ss.utbStrStimAlign});

% pStim trials (pseudoStim trials taken before reach onset)
pStimTrI = cellfun(@(a) ~isempty(a), {ss.utbCgPstimAlign}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCgPstimAlign}, {ss.utbStrPstimAlign});

% noStim trials (with reach and pull)
noStimTrI = cellfun(@(a) isempty(a), {ss.spkTimeBlaserI}) & ...
    cellfun(@(a, b)  excludeNaNsInBothPops(a, b), {ss.unitTimeBCg}, {ss.unitTimeBStr});

valNoStimTrI = noStimTrI & valRchTrI & valPullTrI;
valNoStimTrId = find(valNoStimTrI);

%% RRR with trials without silencing
% build X (n x p matrix containing the residual activity of the source)
% population: M1) and Y (n x q matrix containing the residual activity of the target population: STR) matrices.
%[X.concat, X.numbUnit, X.numbTime, X.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBCtx});
%[Y.concat, Y.numbUnit, Y.numbTime, Y.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBStr});

% X and Y no stim trials aligned to reach start
Xc = {ss(valNoStimTrId).unitTimeBCg};
Yc = {ss(valNoStimTrId).unitTimeBStr};

% X and Y stim trials aligned to laser onset
Xc_stim = {ss(stimTrI).utbCgStimAlign};
Yc_stim = {ss(stimTrI).utbStrStimAlign};

% X and Y psuedo-stim trials aligned to pseudo-laser onset
Xc_pStim = {ss(pStimTrI).utbCgPstimAlign};
Yc_pStim = {ss(pStimTrI).utbStrPstimAlign};

rrrRez_dir = dir(fullfile(saveDir, ['rrrRezCV_stimPstim_Cg_' saveName, sprintf('_Dims%d', dims), sprintf('_Folds%d', folds), '.mat']));
if ~isempty(rrrRez_dir)
    load(fullfile(rrrRez_dir.folder, rrrRez_dir.name), 'rrrCv')
else
    % run cross-validated (trial-shuffled) rrr on reach-aligned data with input-specified dims and folds
    rrrCv = reducedRankRegressCrossVal(Xc, Yc, dims, folds, true);
end

% Test the RRR model on stim and pseudoStim aligned data
for f = 1:size(rrrCv, 2)
    % stim trial test set (cortex silencing trials)
    [rrrCv(f).Xcc_stim, ~, rrrCv(f).numbTimeStim, rrrCv(f).numbTrial_stim] = concatUnitTimeBCell(Xc_stim);
    rrrCv(f).Ycc_stim = concatUnitTimeBCell(Yc_stim);

    % pseudo-stim trial teset set
    [rrrCv(f).Xcc_pStim, ~, rrrCv(f).numbTimePstim, rrrCv(f).numbTrial_pStim] = concatUnitTimeBCell(Xc_pStim);
    rrrCv(f).Ycc_pStim = concatUnitTimeBCell(Yc_pStim);

    % get Yhat_stim (stim trials)
    rrrCv(f).Yhat_stim = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_stim, rrrCv(f).Ycc_stim, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_stim = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_stim, rrrCv(f).numbUnitY, rrrCv(f).numbTimeStim, rrrCv(f).numbTrial_stim);

    % get Yhat_pStim (pseudo-stim trials)
    rrrCv(f).Yhat_pStim = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_pStim, rrrCv(f).Ycc_pStim, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_pStim = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_pStim, rrrCv(f).numbUnitY, rrrCv(f).numbTimePstim, rrrCv(f).numbTrial_pStim);

    % calculate R2
    rrrCv(f).r2_stim = calculateR2(rrrCv(f).Ycc_stim, rrrCv(f).Yhat_stim);
    rrrCv(f).r2_pStim = calculateR2(rrrCv(f).Ycc_pStim, rrrCv(f).Yhat_pStim);
    % calculate R2 trial by trial
    rrrCv(f).r2_stim_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_stim, rrrCv(f).YhatC_stim', 'UniformOutput', false));
    rrrCv(f).r2_pStim_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_pStim, rrrCv(f).YhatC_pStim', 'UniformOutput', false));
end

% organize R2 results
rrrRezR2.mR2_tbyt = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytStim = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_stim_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytPstim = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_pStim_tbyt}, 'UniformOutput', false)));

rrrRezR2.mR2 = nanmean([rrrCv.r2]);
rrrRezR2.mR2_stim = nanmean([rrrCv.r2_stim]);
rrrRezR2.mR2_pStim = nanmean([rrrCv.r2_pStim]);

rrrRezR2.medR2_tbyt = nanmedian(cell2mat(cellfun(@nanmean, {rrrCv.r2_tbyt}, 'UniformOutput', false)));
rrrRezR2.medR2_tbytStim = nanmedian(cell2mat(cellfun(@nanmean, {rrrCv.r2_stim_tbyt}, 'UniformOutput', false)));
rrrRezR2.medR2_tbytPstim = nanmedian(cell2mat(cellfun(@nanmean, {rrrCv.r2_pStim_tbyt}, 'UniformOutput', false)));

rrrRezR2.medR2 = nanmedian([rrrCv.r2]);
rrrRezR2.medR2_stim = nanmedian([rrrCv.r2_stim]);
rrrRezR2.medR2_pStim = nanmedian([rrrCv.r2_pStim]);

%% save results
save(fullfile(saveDir, ['rrrRezCV_stimPstimPrepExtWoTo_Cg_' saveName, sprintf('_Dims%d', dims), sprintf('_Folds%d', folds)]), 'rrrCv', 'rrrRezR2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [concatMat, numbUnit, numbTime, numbTrial] = concatUnitTimeBCell(unitTimeBCell)
        numbTrial = length(unitTimeBCell);

        sizeCell = cellfun(@size, unitTimeBCell, 'UniformOutput', false);
        sizeMat = cell2mat(sizeCell(:));
        assert(length(unique(sizeMat(:, 1)))==1); % ensure that the number of units match
        assert(length(unique(sizeMat(:, 2)))==1); % ensure that the number of time points match

        numbUnit = unique(sizeMat(:, 1));
        numbTime = unique(sizeMat(:, 2));

        concatMat = full(cell2mat(cellfun(@(a) a', unitTimeBCell, 'UniformOutput', false)'));
    end

    function unitTimeBCell = reshapeYhatToUnitTimeBCell(Yhat_concat, numbUnit, numbTime, numbTrial)
        % Orient the input matrix as numbUnit x numbDataPoints (time x trial)
        % because that orientation is required for proper reshaping.
        if size(Yhat_concat, 1) == numbUnit
            Yhat = Yhat_concat;
        elseif size(Yhat_concat, 2) == numbUnit
            Yhat = Yhat_concat';
        end

        % reshape to a 3D array
        rsArray = reshape(Yhat, [numbUnit, numbTime, numbTrial]);

        % convert to cell array
        unitTimeBCell = squeeze(mat2cell(rsArray, numbUnit, numbTime, ones(1, numbTrial)));


    end

    function R2 = calculateR2(Y, Y_hat)
        % Y is the matrix of actual values
        % Y_hat is the matrix of predicted values

        % Ensure Y and Y_hat are the same size
        if size(Y) ~= size(Y_hat)
            error('Y and Y_hat must be the same size');
        end

        % Calculate the total sum of squares (SST)
        SST = sum((Y - mean(Y, 'all')).^2, 'all');

        % Calculate the residual sum of squares (SSR)
        SSR = sum((Y - Y_hat).^2, 'all');

        % Calculate R^2
        R2 = 1 - (SSR / SST);
    end

    function YhatC = getYhatStackedB(X, stackedB)
        %this function computes Yhat by multiplying X with weight matrices stacked
        % over along the 3rd dimension. The stacked B is assumed to have
        % dimensions: # of source units by # of target units by # of stacked weight
        % matrices

        % YhatM
        YhatC = cell(1, 1, size(stackedB, 3));

        % get intercept
        if size(stackedB, 1)-size(X, 2)==0
            intercept = [];
        elseif size(stackedB, 1)-size(X, 2)==1
            intercept = ones(size(X, 1), 1);
        else
            error("Input matrix dimensions do not make sense!")
        end


        for jj = 1:size(stackedB, 3)
            yhat = [intercept X]*stackedB(:, :, jj);
            YhatC{1, 1, jj} = yhat;
        end

    end

    function rrrRez = reducedRankRegressCrossVal(Xc, Yc, dim, folds, trialShuffleLogic)

        % Xc: predictor (a cell array each cell containing neuron by time matrix, e.g., {ss(valNoStimTrId).unitTimeBCtx}).
        % Yc: target (a cell array each cell containing neuron by time matrix, e.g., {ss(valNoStimTrId).unitTimeBStr}).
        % dim: dimension for reduced rank regression
        % folds:
        % trialShuffleLogic:


        assert(length(Xc)==length(Yc))

        [trainIndices, testIndices] = crossValidationFolds(length(Xc), folds, trialShuffleLogic);

        for ff = 1:folds
            % train set
            [rrrRez(ff).Xcc_train, rrrRez(ff).numbUnitX, rrrRez(ff).numbTime, rrrRez(ff).numbTrial_train] = concatUnitTimeBCell(Xc(trainIndices{ff}));
            [rrrRez(ff).Ycc_train, rrrRez(ff).numbUnitY] = concatUnitTimeBCell(Yc(trainIndices{ff}));

            % test set
            [rrrRez(ff).Xcc_test, ~, ~, rrrRez(ff).numbTrial_test] = concatUnitTimeBCell(Xc(testIndices{ff}));
            rrrRez(ff).Ycc_test = concatUnitTimeBCell(Yc(testIndices{ff}));

            % train phase (run RRR)
            [rrrRez(ff).B, rrrRez(ff).B_, rrrRez(ff).V] = ReducedRankRegress(rrrRez(ff).Ycc_train, rrrRez(ff).Xcc_train, ...
                dim, 'RIDGEINIT', true, 'SCALE', true);

            % get Yhat
            rrrRez(ff).Yhat_test = cell2mat(getYhatStackedB(rrrRez(ff).Xcc_test, rrrRez(ff).B)); % yhat for reduced rank
            rrrRez(ff).YhatC_test = reshapeYhatToUnitTimeBCell(rrrRez(ff).Yhat_test, rrrRez(ff).numbUnitY, rrrRez(ff).numbTime, rrrRez(ff).numbTrial_test);

            % calculate R2
            rrrRez(ff).r2 = calculateR2(rrrRez(ff).Ycc_test, rrrRez(ff).Yhat_test);
            % calculate R2 trial by trial
            rrrRez(ff).r2_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc(testIndices{ff}), rrrRez(ff).YhatC_test', 'UniformOutput', false));

            fprintf("Completed rrr cv fold #%d\n", ff)
        end

    end

    function YhatC = getYhatStackedBwithInterceptUpdate(X, Y, stackedB)
        %this function computes Yhat by multiplying X with weight matrices stacked
        % over along the 3rd dimension. The stacked B is assumed to have
        % dimensions: # of source units by # of target units by # of stacked weight
        % matrices

        % YhatM
        YhatC = cell(1, 1, size(stackedB, 3));

        m = mean(X,1);

        % get intercept
        if size(stackedB, 1)-size(X, 2)==0
            intercept = [];
        elseif size(stackedB, 1)-size(X, 2)==1
            intercept = ones(size(X, 1), 1);
        else
            error("Input matrix dimensions do not make sense!")
        end


        for jj = 1:size(stackedB, 3)
            if size(stackedB, 1)-size(X,2)==1
                B = stackedB(:, :, jj);
                % Note that the intercept term needs to be updated especially when
                % applying the weight matrix to a different epoch, where the mean
                % differs.
                B_correct = [ mean(Y,1) - m*B(2:end, :); B(2:end, :)];
                yhat = [intercept X]*B_correct;
                YhatC{1, 1, jj} = yhat;
            end

        end

    end
end



