function accuracy = multiClass_naiveBayesClassifier_peth(Xs, y, resample)
% Xs = [cell2mat(hitFstLickDffSsItp); cell2mat(faFstLickDffSsItp); cell2mat(postStimLickDffSsItp); cell2mat(itiLickDffSsItp)];
% y = [ones(size(hitFstLickDffSsItp, 1), 1)*1; ones(size(faFstLickDffSsItp, 1), 1)*2; ones(size(postStimLickDffSsItp, 1), 1)*3; ones(size(itiLickDffSsItp, 1), 1)*4];
% resample = 10; 


accuracy = zeros(resample, size(Xs, 2)); 

% Find the minimum class size
minSize = min(histcounts(y));

for tt = 1:size(Xs, 2)

    X = Xs(:, tt);

    for rs = 1:resample
        % Create new variables for balanced data
        X_balanced = [];
        y_balanced = [];

        % Balance each class
        uniqueClasses = unique(y);
        for i = 1:length(uniqueClasses)
            classIndices = find(y == uniqueClasses(i));
            % Randomly select 'minSize' samples from each class
            randIndices = randsample(classIndices, minSize);
            X_balanced = [X_balanced; X(randIndices, :)];
            y_balanced = [y_balanced; y(randIndices)];
        end

        % Split Data into Training and Testing
        cv = cvpartition(size(X_balanced,1),'HoldOut',0.3);
        idx = cv.test;

        % Separate to training and test data
        XTrain = X_balanced(~idx,:);
        YTrain = y_balanced(~idx,:);
        XTest  = X_balanced(idx,:);
        YTest  = y_balanced(idx,:);

        % Step 4: Train SVM Classifier
        %nbModel = fitcecoc(XTrain, YTrain);
        nbModel = fitcnb(XTrain, YTrain);

        % Step 5: Test the Classifier
        YPred = predict(nbModel, XTest);

        % Evaluate performance
        accuracy(rs, tt) = sum(YPred == YTest) / length(YTest);
        fprintf('Finished resample #%d of bin #%d with accuracy: %.2f%%\n', rs, tt, accuracy(rs, tt) * 100);
    end
end

end

