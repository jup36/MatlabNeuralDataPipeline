function [rrrCv, rrrRezR2] = stimEffect_neuralTrj_cv_stim_pstim_rch_noRch(filePath, dims, folds)
%This function is to run the reduced rank regression (RRR) with cross-validation,
% which trains the RRR model that predicts the target (e.g., str) neural population activity
% based on the source neural population activity (e.g., ctx) with cross validation.
% The trained models are tested both on the held-out no-stim and stim
% trials. The results of rrr and testing (r2) are saved and/or returned as
% an output of this function.

saveNameIdx = strfind(filePath, 'WR');
saveName = filePath(saveNameIdx(1):saveNameIdx(1)+10);
saveDir = fileparts(filePath);

%filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/js2p0_tbytSpkHandJsTrjBin_50ms_WR40_081919.mat';
load(fullfile(filePath), 'ss')

%% get indices
% valid reach and pull trials (trials with valid hand trajectories)
valRchTrI = cellfun(@(a) ~isempty(a), {ss.spkRchIdx});
valPullTrI = cellfun(@(a) ~isempty(a), {ss.spkPullIdx});

excludeNaNsInBothPops = @(a, b) sum(sum(isnan(full(a)), 2))==0 & sum(sum(isnan(full(b)), 2))==0;

% stim reach trials
stimRchTrI =  cellfun(@(a) ~isempty(a), {ss.utbCtxStimRch}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCtxStimRch}, {ss.utbStrStimRch});

% pStim reach trials (pseudoStim trials taken before reach onset aligned to the timing at which stim would've delivered)
pStimRchTrI = cellfun(@(a) ~isempty(a), {ss.utbCtxPstimRch}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCtxPstimRch}, {ss.utbStrPstimRch});

% stim no-reach trials
stimNoRchTrI = cellfun(@(a) ~isempty(a), {ss.utbCtxStimNoRch}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCtxStimNoRch}, {ss.utbStrStimNoRch});

% stim no-reach trials
pStimNoRchTrI = cellfun(@(a) ~isempty(a), {ss.utbCtxPstimNoRch}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCtxPstimNoRch}, {ss.utbStrPstimNoRch});


% noStim trials (with reach and pull)
noStimTrI = cellfun(@(a) isempty(a), {ss.spkTimeBlaserI}) & ...
    cellfun(@(a, b)  excludeNaNsInBothPops(a, b), {ss.unitTimeBCtx}, {ss.unitTimeBStr});

valNoStimTrI = noStimTrI & valRchTrI & valPullTrI;
valNoStimTrId = find(valNoStimTrI);

%% RRR with trials without silencing
% build X (n x p matrix containing the residual activity of the source)
% population: M1) and Y (n x q matrix containing the residual activity of the target population: STR) matrices.
%[X.concat, X.numbUnit, X.numbTime, X.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBCtx});
%[Y.concat, Y.numbUnit, Y.numbTime, Y.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBStr});

% X and Y no stim trials aligned to reach start
Xc = {ss(valNoStimTrId).unitTimeBCtx};
Yc = {ss(valNoStimTrId).unitTimeBStr};

% X and Y stim reach trials aligned to reach onset
Xc_stimRch = {ss(stimRchTrI).utbCtxStimRch};
Yc_stimRch = {ss(stimRchTrI).utbStrStimRch};

% X and Y stim NO reach trials aligned to stim onset
Xc_stimNoRch = {ss(stimNoRchTrI).utbCtxStimNoRch};
Yc_stimNoRch = {ss(stimNoRchTrI).utbStrStimNoRch};

% X and Y pStim reach trials aligned to reach onset
Xc_pStimRch = {ss(pStimRchTrI).utbCtxPstimRch};
Yc_pStimRch = {ss(pStimRchTrI).utbStrPstimRch};

% X and Y pStim NO reach trials aligned to pStim onset
Xc_pStimNoRch = {ss(pStimNoRchTrI).utbCtxPstimNoRch};
Yc_pStimNoRch = {ss(pStimNoRchTrI).utbStrPstimNoRch};

rrrRez_dir = dir(fullfile(saveDir, ['rrrRezCV_stimPstim_' saveName, sprintf('_Dims%d', dims), sprintf('_Folds%d', folds), '.mat']));
if ~isempty(rrrRez_dir)
    load(fullfile(rrrRez_dir.folder, rrrRez_dir.name), 'rrrCv')
else
    % run cross-validated (trial-shuffled) rrr on reach-aligned data with input-specified dims and folds
    rrrCv = reducedRankRegressCrossVal(Xc, Yc, dims, folds, true);
end

% Test the RRR model on stim and pseudoStim aligned data
for f = 1:size(rrrCv, 2)
    % stim reach trial test set
    [rrrCv(f).Xcc_stimRch, ~, rrrCv(f).numbTimeStimRch, rrrCv(f).numbTrial_stimRch] = concatUnitTimeBCell(Xc_stimRch);
    rrrCv(f).Ycc_stimRch = concatUnitTimeBCell(Yc_stimRch);

    % stim no reach trial test set
    [rrrCv(f).Xcc_stimNoRch, ~, rrrCv(f).numbTimeStimNoRch, rrrCv(f).numbTrial_stimNoRch] = concatUnitTimeBCell(Xc_stimNoRch);
    rrrCv(f).Ycc_stimNoRch = concatUnitTimeBCell(Yc_stimNoRch);

    % pStim reach trial test set
    [rrrCv(f).Xcc_pStimRch, ~, rrrCv(f).numbTimePstimRch, rrrCv(f).numbTrial_pStimRch] = concatUnitTimeBCell(Xc_pStimRch);
    rrrCv(f).Ycc_pStimRch = concatUnitTimeBCell(Yc_pStimRch);

    % pStim no reach trial test set
    [rrrCv(f).Xcc_pStimNoRch, ~, rrrCv(f).numbTimePstimNoRch, rrrCv(f).numbTrial_pStimNoRch] = concatUnitTimeBCell(Xc_pStimNoRch);
    rrrCv(f).Ycc_pStimNoRch = concatUnitTimeBCell(Yc_pStimNoRch);

    % get Yhat_stimRch (stimRch trials)
    rrrCv(f).Yhat_stimRch = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_stimRch, rrrCv(f).Ycc_stimRch, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_stimRch = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_stimRch, rrrCv(f).numbUnitY, rrrCv(f).numbTimeStimRch, rrrCv(f).numbTrial_stimRch);

    % get Yhat_stimNoRch (stimNoRch trials)
    rrrCv(f).Yhat_stimNoRch = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_stimNoRch, rrrCv(f).Ycc_stimNoRch, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_stimNoRch = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_stimNoRch, rrrCv(f).numbUnitY, rrrCv(f).numbTimeStimNoRch, rrrCv(f).numbTrial_stimNoRch);

    % get Yhat_pStimRch (pStimRch trials)
    rrrCv(f).Yhat_pStimRch = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_pStimRch, rrrCv(f).Ycc_pStimRch, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_pStimRch = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_pStimRch, rrrCv(f).numbUnitY, rrrCv(f).numbTimePstimRch, rrrCv(f).numbTrial_pStimRch);

    % get Yhat_pStimNoRch (pStimNoRch trials)
    rrrCv(f).Yhat_pStimNoRch = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_pStimNoRch, rrrCv(f).Ycc_pStimNoRch, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_pStimNoRch = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_pStimNoRch, rrrCv(f).numbUnitY, rrrCv(f).numbTimePstimNoRch, rrrCv(f).numbTrial_pStimNoRch);

    % calculate R2
    rrrCv(f).r2_stimRch = calculateR2(rrrCv(f).Ycc_stimRch, rrrCv(f).Yhat_stimRch);
    rrrCv(f).r2_stimNoRch = calculateR2(rrrCv(f).Ycc_stimNoRch, rrrCv(f).Yhat_stimNoRch);
    rrrCv(f).r2_pStimRch = calculateR2(rrrCv(f).Ycc_pStimRch, rrrCv(f).Yhat_pStimRch);
    rrrCv(f).r2_pStimNoRch = calculateR2(rrrCv(f).Ycc_pStimNoRch, rrrCv(f).Yhat_pStimNoRch);

    % calculate R2 trial by trial
    rrrCv(f).r2_stimRch_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_stimRch, rrrCv(f).YhatC_stimRch', 'UniformOutput', false));
    rrrCv(f).r2_stimNoRch_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_stimNoRch, rrrCv(f).YhatC_stimNoRch', 'UniformOutput', false));
    rrrCv(f).r2_pStimRch_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_pStimRch, rrrCv(f).YhatC_pStimRch', 'UniformOutput', false));
    rrrCv(f).r2_pStimNoRch_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_pStimNoRch, rrrCv(f).YhatC_pStimNoRch', 'UniformOutput', false));
end

% organize R2 results
rrrRezR2.mR2_tbyt = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytStimRch = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_stimRch_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytStimNoRch = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_stimNoRch_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytPstimRch = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_pStimRch_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytPstimNoRch = nanmean(cell2mat(cellfun(@nanmean, {rrrCv.r2_pStimNoRch_tbyt}, 'UniformOutput', false)));

rrrRezR2.mR2 = nanmean([rrrCv.r2]);
rrrRezR2.mR2_stimRch = nanmean([rrrCv.r2_stimRch]);
rrrRezR2.mR2_stimNoRch = nanmean([rrrCv.r2_stimNoRch]);
rrrRezR2.mR2_pStimRch = nanmean([rrrCv.r2_pStimRch]);
rrrRezR2.mR2_pStimNoRch = nanmean([rrrCv.r2_pStimNoRch]);

rrrRezR2.medR2_tbyt = nanmedian(cell2mat(cellfun(@nanmedian, {rrrCv.r2_tbyt}, 'UniformOutput', false)));
rrrRezR2.medR2_tbytStimRch = nanmedian(cell2mat(cellfun(@nanmedian, {rrrCv.r2_stimRch_tbyt}, 'UniformOutput', false)));
rrrRezR2.medR2_tbytStimNoRch = nanmedian(cell2mat(cellfun(@nanmedian, {rrrCv.r2_stimNoRch_tbyt}, 'UniformOutput', false)));
rrrRezR2.medR2_tbytPstimRch = nanmedian(cell2mat(cellfun(@nanmedian, {rrrCv.r2_pStimRch_tbyt}, 'UniformOutput', false)));
rrrRezR2.medR2_tbytPstimNoRch = nanmedian(cell2mat(cellfun(@nanmedian, {rrrCv.r2_pStimNoRch_tbyt}, 'UniformOutput', false)));

rrrRezR2.medR2 = nanmedian([rrrCv.r2]);
rrrRezR2.medR2_stimRch = nanmedian([rrrCv.r2_stimRch]);
rrrRezR2.medR2_stimNoRch = nanmedian([rrrCv.r2_stimNoRch]);
rrrRezR2.medR2_pStimRch = nanmedian([rrrCv.r2_pStimRch]);
rrrRezR2.medR2_pStimNoRch = nanmedian([rrrCv.r2_pStimNoRch]);

%% save results
save(fullfile(saveDir, ['rrrRezCV_stimPstimPrep_' saveName, sprintf('_Dims%d', dims), sprintf('_Folds%d', folds)]), 'rrrCv', 'rrrRezR2')

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



