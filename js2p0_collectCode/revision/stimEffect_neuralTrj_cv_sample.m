function [rrrCv, rrrRezR2] = stimEffect_neuralTrj_cv_sample(filePath, dims, folds)
%This function is to run the reduced rank regression (RRR) with cross-validation,
% which trains the RRR model that predicts the target (e.g., str) neural population activity
% based on the source neural population activity (e.g., ctx) with cross validation.
% The trained models are tested both on the held-out no-stim and stim
% trials. The results of rrr and testing (r2) are saved and/or returned as
% an output of this function.

%filePath = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/js2p0_tbytSpkHandJsTrjBin_50ms_WR40_082019.mat';
load(fullfile(filePath), 'ss')

%% get indices
% valid reach and pull trials (trials with valid hand trajectories)
valRchTrI = cellfun(@(a) ~isempty(a), {ss.spkRchIdx});
valPullTrI = cellfun(@(a) ~isempty(a), {ss.spkPullIdx});

% stim trials
stimTrI = cellfun(@(a) ~isempty(a), {ss.spkTimeBlaserI});
valStimTrI = stimTrI & valRchTrI & valPullTrI;
valStimTrId = find(valStimTrI);

% noStim trials
noStimTrI = cellfun(@(a) isempty(a), {ss.spkTimeBlaserI});
valNoStimTrI = noStimTrI & valRchTrI & valPullTrI;
valNoStimTrId = find(valNoStimTrI);

stbLaserIC = {ss(valStimTrI).spkTimeBlaserI}; % spikeTimeBin laser index
stbRchIC = {ss(valStimTrI).spkRchIdx}; % spikeTimeBin reach index
stbPullIC = {ss(valStimTrI).spkPullIdx}; % spikeTimeBin pull index

% check overlap between stim and reach and pull
%laserRchIC = cellfun(@(a, b) a & b, stbLaserIC, stbRchIC, 'UniformOutput', false);
%cellfun(@sum, laserRchIC)

%laserPullIC = cellfun(@(a, b) a & b, stbLaserIC, stbPullIC, 'UniformOutput', false);
%cellfun(@sum, laserPullIC)

%% RRR with trials without silencing
% build X (n x p matrix containing the residual activity of the source)
% population: M1) and Y (n x q matrix containing the residual activity of the target population: STR) matrices.
[X.concat, X.numbUnit, X.numbTime, X.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBCtx});
[Y.concat, Y.numbUnit, Y.numbTime, Y.numbTrial] = concatUnitTimeBCell({ss(valNoStimTrId).unitTimeBStr});

Xc = {ss(valNoStimTrId).unitTimeBCtx};
Yc = {ss(valNoStimTrId).unitTimeBStr};

Xc_stim = {ss(valStimTrId).unitTimeBCtx};
Yc_stim = {ss(valStimTrId).unitTimeBStr};

% run cross-validated (trial-shuffled) rrr with 10 dimensions and 10 folds
rrrCv = reducedRankRegressCrossVal(Xc, Yc, dims, folds, true);

% Test the RRR model on stim trial data
for f = 1:size(rrrCv, 2)

    numbSample = rrrCv(f).numbTrial_test; 
    stimSamples = randperm(length(Xc_stim), numbSample); 

    % test set (cortex silencing trials)
    [rrrCv(f).Xcc_stim, ~, ~, rrrCv(f).numbTrial_stim] = concatUnitTimeBCell(Xc_stim(stimSamples));
    rrrCv(f).Ycc_stim = concatUnitTimeBCell(Yc_stim(stimSamples));

    % get Yhat_stim
    rrrCv(f).Yhat_stim = cell2mat(getYhatStackedB(rrrCv(f).Xcc_stim, rrrCv(f).B)); % yhat for reduced rank
    rrrCv(f).YhatC_stim = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_stim, rrrCv(f).numbUnitY, rrrCv(f).numbTime, rrrCv(f).numbTrial_stim);

    % calculate R2
    rrrCv(f).r2_stim = calculateR2(rrrCv(f).Ycc_stim, rrrCv(f).Yhat_stim);
    % calculate R2 trial by trial
    rrrCv(f).r2_stim_tbyt = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_stim(stimSamples), rrrCv(f).YhatC_stim', 'UniformOutput', false));
end

rrrRezR2.mR2_tbyt = mean(cell2mat(cellfun(@mean, {rrrCv.r2_tbyt}, 'UniformOutput', false)));
rrrRezR2.mR2_tbytStim = mean(cell2mat(cellfun(@mean, {rrrCv.r2_stim_tbyt}, 'UniformOutput', false)));

rrrRezR2.mR2 = mean([rrrCv.r2]);
rrrRezR2.mR2_stim = mean([rrrCv.r2_stim]);

%% save results
saveNameIdx = strfind(filePath, 'WR');
saveName = filePath(saveNameIdx(1):saveNameIdx(1)+10);
saveDir = fileparts(filePath);

save(fullfile(saveDir, ['rrrRezCV_' saveName, sprintf('_Dims%d', dims), sprintf('_Folds%d', folds)]), 'rrrCv', 'rrrRezR2')


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

end




