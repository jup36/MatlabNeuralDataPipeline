function PS3_4
%% %%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% cd('C:\Users\jup36\Dropbox\NSP\PS3')      % the corresponding directory
load('ps3_realdata.mat','-mat');
[NumTrainData NumClass]=size(train_trial);
NumTestData=size(test_trial,1);
for classIX=1:NumClass
    for trainDataIX=1:NumTrainData
        currData=train_trial(trainDataIX,classIX).spikes(:,351:550);
        trainDataArr(classIX,trainDataIX,:)=sum(currData,2);
    end
    for testDataIX=1:NumTestData
        currData=test_trial(testDataIX,classIX).spikes(:,351:550);
        testDataArr(classIX,testDataIX,:)=sum(currData,2);
    end
end
NumFea=size(trainDataArr,3);
% For test data
actLabel=repmat([1:NumClass]',1,NumTestData);
testData=reshape(testDataArr,[],NumFea);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Part (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% Fit the parameters of model(1)
modMean=squeeze(mean(trainDataArr,2));
% Remove the mean of each class
trainDataArrRemMean=reshape(repmat(modMean,1,NumTrainData),[NumClass NumFea NumTrainData]);
trainDataArrRemMean=trainDataArr-permute(trainDataArrRemMean,[1 3 2]);
% Shared covariance matrix
modCov{1}=cov(reshape(trainDataArrRemMean,[],NumFea));
% Test
for classIX=1:NumClass...
    tmp=testData-repmat(modMean(classIX,:),NumTestData*NumClass,1);
    logP(:,classIX)=sum(tmp*inv(modCov{1}).*tmp,2);
end
[minVal predLabel]=min(logP,[],2);
corrPred=find((predLabel-actLabel(:))==0);
corrRatio_a=length(corrPred)/(NumTestData*NumClass);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Part (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% Fit the parameters of model(2)
% The covariance matrices are singular
modMean=squeeze(mean(trainDataArr,2));
for classIX=1:NumClass
    modCov{classIX}=cov(squeeze(trainDataArr(classIX,:,:)));
end
% Test
for classIX=1:NumClass
    tmp=testData-repmat(modMean(classIX,:),NumTestData*NumClass,1);
    logP(:,classIX)=sum(tmp*inv(modCov{classIX}).*tmp,2)+log(det(modCov{classIX}));
end
[minVal predLabel]=min(logP,[],2);corrPred=find((predLabel-actLabel(:))==0);
corrRatio_b=length(corrPred)/(NumTestData*NumClass);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Part (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% Fit the parameters of model (3)
modMean=squeeze(mean(trainDataArr,2));
% Test
for classIX=1:NumClass
    currLambda=modMean(classIX,:);
    logP(:,classIX)=-testData*log(currLambda')+sum(currLambda);
end
[minVal predLabel]=min(logP,[],2);
corrPred=find((predLabel-actLabel(:))==0);
corrRatio_c=length(corrPred)/(NumTestData*NumClass);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Part (d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% Fit the parameters of model (3) with min variance set
modMean=squeeze(mean(trainDataArr,2));
regCoef=0.01;
modMean=max(modMean,regCoef);
% Test
for classIX=1:NumClass
    currLambda=modMean(classIX,:);
    logP(:,classIX)=-testData*log(currLambda')+sum(currLambda);
end
[minVal predLabel]=min(logP,[],2);
corrPred=find((predLabel-actLabel(:))==0);
corrRatio_d=length(corrPred)/(NumTestData*NumClass);
%% %%%%%%%%%%%%% Print out the classification results %%%%%%%%%%%%%%%%%
%%%
fprintf('(a) Classification accuracy of model(1): %2.2f%
%\n',corrRatio_a*100);
fprintf('(b) Classification accuracy of model(2): %2.2f%
%\n',corrRatio_b*100);
fprintf(' Not enough training data to fit a full covariance matrix
(random chance level)\n');
fprintf('(c) Classification accuracy of model(3): %2.2f%
%\n',corrRatio_c*100);
fprintf(' Poisson distribution has 0 var for \mu=0 (random chance
level)\n');
fprintf('(d) Classification accuracy of model(3) with minimum variance
set: %2.2f%%\n',corrRatio_d*100);
return;
(a) Classification accuracy of model(1): 96.02%
(b) Classification accuracy of model(2): 12.5%
Not enough training data to fit a full covariance matrix (random
chance level)
(c) Classification accuracy of model(3): 12.5%
Poisson distribution has 0 var for \mu=0 (random chance level)
(d) Classification accuracy of model(3) with minimum variance set:
94.09%