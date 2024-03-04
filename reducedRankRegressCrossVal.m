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
    
    numB = size(rrrRez(ff).B, 2)/rrrRez(ff).numbUnitY; % the number of weight matrix B; it must equal to the number of dims
    rrrRez(ff).stackedB = zeros(size(rrrRez(ff).B, 1), rrrRez(ff).numbUnitY, numB); 
    if numB==length(dim)
        if numB >=2
            for dd = 0:numB-1
                rrrRez(ff).stackedB(:, :, dd+1) = rrrRez(ff).B(:, dd*rrrRez(ff).numbUnitY+1:(dd+1)*rrrRez(ff).numbUnitY); 
            end
        elseif numB == 1
            rrrRez(ff).stackedB = rrrRez(ff).B; 
        end
    else
        error("The dimension of weight matrix B doesn't make sense, check dimensions!")
    end

    % get Yhat
    rrrRez(ff).Yhat_test = cell2mat(getYhatStackedB(rrrRez(ff).Xcc_test, rrrRez(ff).stackedB)); % yhat for reduced rank 
 
    for jj = 1:size(rrrRez(ff).Yhat_test, 3) % increment dims
        % calculate R2 over all trials: rrrRez(ff).r2 will be dims-by-1 
        rrrRez(ff).r2(jj, 1) = calculateR2(rrrRez(ff).Ycc_test, rrrRez(ff).Yhat_test(:, :, jj)); 
        % calculate R2 over each trial: rrrRez(ff).r2_tbyt will be dims-by-trials
        rsYhatC_test = reshapeYhatToUnitTimeBCell(rrrRez(ff).Yhat_test(:, :, jj), rrrRez(ff).numbUnitY, rrrRez(ff).numbTime, rrrRez(ff).numbTrial_test); 
        rrrRez(ff).r2_tbyt(jj, :) = cell2mat(cellfun(@(a, b) calculateR2(a, b), Yc(testIndices{ff}), rsYhatC_test', 'UniformOutput', false)); 
    end
    fprintf("Completed rrr cv fold #%d\n", ff)
end

end

