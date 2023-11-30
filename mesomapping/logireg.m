function accuracy = logireg(X, Y)
% % X: observations-by-feature
% %   e.g., cell2mat(stimOnDff_m1_itp); 
% % Y: label 

Y = Y(:);

% Determine the number of samples for each class in the test set
numClass0 = sum(Y == 0);
numClass1 = sum(Y == 1);
testSetSize = min(numClass0, numClass1) / 2; % Or any other criterion

% Randomly sample for the test set
indices0 = find(Y == 0);
indices1 = find(Y == 1);
testIndices0 = randsample(indices0, testSetSize);
testIndices1 = randsample(indices1, testSetSize);

% Combine test samples and create the test set
testInd = [testIndices0; testIndices1];
testInd = testInd(randperm(length(testInd)));

% Use the remaining data as the training set
trainInd = setdiff(1:length(Y), testInd);

accuracy = nan(1, size(X, 2)); 
% Training
for jj = 1:size(X, 2) % iterate features
    model = fitglm(X(trainInd, jj), Y(trainInd), 'Distribution', 'binomial');
    % Testing
    Y_pred = predict(model, X(testInd,jj)) > 0.5; % Thresholding at 0.5
    accuracy(1, jj) = sum(Y_pred == Y(testInd)) / length(Y_pred);
end

end

