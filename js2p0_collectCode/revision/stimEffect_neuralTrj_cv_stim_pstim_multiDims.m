function [rrrCv, rrrRezR2] = stimEffect_neuralTrj_cv_stim_pstim_multiDims(filePath, folds)
%This function is to run the reduced rank regression (RRR) with cross-validation,
% which trains the RRR model that predicts the target (e.g., str) neural population activity
% based on the source neural population activity (e.g., ctx) with cross validation.
% The trained models are tested both on the held-out no-stim and stim
% trials. The results of rrr and testing (r2) are saved and/or returned as
% an output of this function.
%filePath = '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/js2p0_tbytSpkHandJsTrjBin_50ms_stimPstim_WR40_081919.mat';

saveNameIdx = strfind(filePath, 'WR');
saveName = filePath(saveNameIdx(1):saveNameIdx(1)+10);
saveDir = fileparts(filePath);

load(fullfile(filePath), 'ss')

%% get indices
% valid reach and pull trials (trials with valid hand trajectories)
valRchTrI = cellfun(@(a) ~isempty(a), {ss.spkRchIdx});
valPullTrI = cellfun(@(a) ~isempty(a), {ss.spkPullIdx});

excludeNaNsInBothPops = @(a, b) sum(sum(isnan(full(a)), 2))==0 & sum(sum(isnan(full(b)), 2))==0;


% stim trials
stimTrI =  cellfun(@(a) ~isempty(a), {ss.utbCtxStimAlign}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCtxStimAlign}, {ss.utbStrStimAlign});

% pStim trials (pseudoStim trials taken before reach onset)
pStimTrI = cellfun(@(a) ~isempty(a), {ss.utbCtxPstimAlign}) & ...
    cellfun(@(a, b) excludeNaNsInBothPops(a, b), {ss.utbCtxPstimAlign}, {ss.utbStrPstimAlign});

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

% X and Y stim trials aligned to laser onset
Xc_stim = {ss(stimTrI).utbCtxStimAlign};
Yc_stim = {ss(stimTrI).utbStrStimAlign};

% X and Y psuedo-stim trials aligned to pseudo-laser onset
Xc_pStim = {ss(pStimTrI).utbCtxPstimAlign};
Yc_pStim = {ss(pStimTrI).utbStrPstimAlign};

dims = [1:20, unique(cellfun(@(a) size(a, 1), Yc))]; % test dimensions (end: full)

rrrRez_dir = dir(fullfile(saveDir, ['rrrRezCV_stimPstim_multiDim_' saveName, sprintf('_Folds%d', folds), '.mat']));
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

    % get Yhat_stim (stim trials), Yhat_pStim (pseudo-stim trials)
    rrrCv(f).Yhat_stim = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_stim, rrrCv(f).Ycc_stim, rrrCv(f).stackedB)); % yhat for reduced rank
    rrrCv(f).Yhat_pStim = cell2mat(getYhatStackedBwithInterceptUpdate(rrrCv(f).Xcc_pStim, rrrCv(f).Ycc_pStim, rrrCv(f).stackedB)); % yhat for reduced rank

    % calculate r2 stim and pStim
    for dd = 1:size(rrrCv(f).Yhat_stim, 3)
        % calculate R2 over all trials: rrrRez(f).r2_* will be dims-by-1
        rrrCv(f).r2_stim(dd, 1) = calculateR2(rrrCv(f).Ycc_stim, rrrCv(f).Yhat_stim(:, :, dd));
        rrrCv(f).r2_pStim(dd, 1) = calculateR2(rrrCv(f).Ycc_pStim, rrrCv(f).Yhat_pStim(:, :, dd));
        % calculate R2 over each trial: rrrRez(ff).r2_*_tbyt will be dims-by-trials
        rsYhatC_stim = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_stim(:, :, dd), rrrCv(f).numbUnitY, rrrCv(f).numbTimeStim, rrrCv(f).numbTrial_stim);
        rrrCv(f).r2_stim_tbyt(dd, :) = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_stim, rsYhatC_stim', 'UniformOutput', false));

        rsYhatC_pStim = reshapeYhatToUnitTimeBCell(rrrCv(f).Yhat_pStim(:, :, dd), rrrCv(f).numbUnitY, rrrCv(f).numbTimePstim, rrrCv(f).numbTrial_pStim);
        rrrCv(f).r2_pStim_tbyt(dd, :) = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc_pStim, rsYhatC_pStim', 'UniformOutput', false));
    end
end

% organize R2 results
rrrRezR2.R2_tbyt = {rrrCv.r2_tbyt};
rrrRezR2.R2_tbytStim = {rrrCv.r2_stim_tbyt};
rrrRezR2.R2_tbytPstim = {rrrCv.r2_pStim_tbyt};

rrrRezR2.R2 = [rrrCv.r2];
rrrRezR2.R2_stim = [rrrCv.r2_stim];
rrrRezR2.R2_pStim = [rrrCv.r2_pStim];

% plot R2 across dimensions relative to the mean full dimensional R2
[rrrRezR2.meanR2, ~, rrrRezR2.semR2] = meanstdsem(rrrRezR2.R2');
meanSemErrorbar(rrrRezR2.meanR2(1:20), rrrRezR2.semR2(1:20)); % plot dim 1 to 20
hold on;
plot(1:20, repmat(mean(rrrRezR2.R2(end, :)), 1, 20), ':') % plot mean full dim result
set(gca, 'TickDir', 'out')

print(fullfile(fileparts(filePath), 'Figure', 'meanSEM_R2_multiDim'), '-vector', '-dpdf')

%% save results
save(fullfile(saveDir, ['rrrRezCV_stimPstim_multiDim_' saveName, sprintf('_Folds%d', folds)]), 'rrrCv', 'rrrRezR2')

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

    function fig = meanSemErrorbar(mean_values, sem_values)

        fig = figure;
        % Assuming mean_values and sem_values are already defined as 1x20 vectors

        % Generate a vector for the x-axis
        x_values = 1:length(mean_values);

        % Plot filled circles for mean values
        scatter(x_values, mean_values, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');

        % Hold on to the current figure
        hold on;

        % Add error bars for standard error of the mean
        errorbar(x_values, mean_values, sem_values, 'LineStyle', 'none', 'Color', 'r', 'CapSize', 10);

        % Label the axes
        xlabel('Dims');
        ylabel('Mean Value');

        % Add a title and a legend
        %title('Mean Values with SEM');
        %legend('Mean values', 'SEM', 'Location', 'best');

        % Hold off to finish the plotting
        hold off;

    end

end



