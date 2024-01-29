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

