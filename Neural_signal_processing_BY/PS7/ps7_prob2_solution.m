
clear all; close all; clear functions; clc

%%
% addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps, modify these!
%filedirectory = 'C:\Users\jup36\Dropbox\NSP\PS7';    % labPC directory  
filedirectory = 'C:\Users\Junchol\Dropbox\NSP\PS7';    % labtop directory  
cd(filedirectory)       
load('ps7_data.mat');

[NTime, NSpikes]=size(Spikes);  % # of time points, # of spikes 

% Project the spikes to a two-dimensional space using PCA
[comp, score]=princomp(Spikes');

% negate the directios
comp = -comp;
score = -score;
NumFea=2;   % # of features (PC1, PC2)
selFea=score(:,1:NumFea);   % selected features: PC1 score, PC2 score
% Compute mean of the spikes
meanSpikes=mean(Spikes,2);

%% GMM model selection
% Set the K= 1, ..., 8
NumClusterList=1:8;
plotCount=0;
for clusterIX=1:length(NumClusterList)
    NumCluster=NumClusterList(clusterIX);
    % Initialize the GMM parameters
    currInitParams.mu=InitParams.mu(:,1:NumCluster);
    for k=1:NumCluster
        currInitParams.Sigma(:,:,k)=InitParams.Sigma;
        currInitParams.pi(k)=1/NumCluster;
    end
    % Corss-validation
    NumFold=4;
    NumPerFold=NSpikes/NumFold;
    for foldIX=1:NumFold
        % Separate the training set and the test set
        testIX=((foldIX-1)*NumPerFold+1):(foldIX*NumPerFold);
        trainIX=setdiff(1:NSpikes,testIX);
        testData=selFea(testIX,:);
        trainData=selFea(trainIX,:);
        % Training
        [mu, sigma, ppi]=func_GMM(currInitParams, trainData');
        % Test (compute the likelihood of each test fold)
        const = -0.5 * NumFea * log(2*pi);
        logMat = nan(NumCluster, length(testIX));
        for k = 1:NumCluster
            S = sigma(:,:,k);
            xdif = bsxfun(@minus, testData', mu(:,k));
            term1 = -0.5 * sum((xdif' * inv(S)) .* xdif', 2); % N x 1
            term2 = const - 0.5 * log(det(S)) + log(ppi(k)); % scalar
            logMat(k,:) = term1' + term2;
        end
        astar = max(logMat, [], 1);
        adif = bsxfun(@minus, logMat, astar);
        nLL = log(sum(exp(adif), 1)) + astar; % 1 x N
        LL(foldIX) = sum(nLL);
        %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (b) %%%%%%%%%%%%%%%%%%%%%%%%
        % Plot the Gaussian distribution from the first CV fold
        if (foldIX==1)
            plotCount=plotCount+1;
            figure(plotCount);
            % Plot the training data in blue
            plot(trainData(:,1),trainData(:,2),'b.');hold on;
            % Plot the training data in red
            plot(testData(:,1),testData(:,2),'r.'); hold on;
            for k=1:NumCluster
                func_plotEllipse(mu(:,k),sigma(:,:,k));
                hold on;
            end
            xlabel('PC1 score'); ylabel('PC2 score');
            titleStr=sprintf('ClusterNum=%d',NumCluster);
            title(titleStr);
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (c) %%%%%%%%%%%%%%%%%%%%%%%%
        % For K=3, plot the waveforms corresponding to GMM cluster centers
        if (NumCluster==3 && foldIX==1)
            recCenter=comp(:,1:2)*mu+repmat(meanSpikes,1,3);
            plotCount=plotCount+1;
            figure(plotCount); plot(recCenter);
            xlabel('Sample'); ylabel('Voltage (\muV)');
            axis tight
            title ('K=3: Waveforms corresponding to GMM cluster centers');
        end
    end
    % The CV likelihoods for K=NumClusterList(clusterIX)
    totLL(clusterIX)=sum(LL);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% Part (a) %%%%%%%%%%%%%%%%%%%%%%%%
% Plot the CV likelihoods versus K
plotCount=plotCount+1;
figure(plotCount); plot(totLL,'-*')
xlabel('K'); ylabel('Cross-validated likelihoods');