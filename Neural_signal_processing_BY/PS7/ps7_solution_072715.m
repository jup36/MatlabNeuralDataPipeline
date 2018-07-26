%% Problem 7.1
clear all; close all; clear functions; clc

%%
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps, modify these!
%filedirectory = 'C:\Users\jup36\Dropbox\NSP\PS7';    % labPC directory  
filedirectory = 'C:\Users\Junchol\Dropbox\NSP\PS7';    % labtop directory  
cd(filedirectory)       
load('ps7_data.mat');
[NTime, NSpikes] = size(Spikes);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the raw spikes
% subplot(2,2,1); plot(Spikes);
% xlabel('Sample'); ylabel('Voltage (\muV)');
% axis tight

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply PCA
% S = cov(Spikes');     % covariance matrix 
% [eigvec, eigval] = eig(S);    % eigen decomposition; coeffs are eigenvalues; scores are eigenvectors; S*v = v*D   

% Project the spikes to a two-dimensional space using PCA
[comp, score]=pca(Spikes');

% negate the directios
comp = -comp;
score = -score;

% Plot the first three principal components
% subplot(2, 2, 2); plot(comp(:, 1), '-r');
% hold on; plot(comp(:, 2), '-g'); hold on; plot(comp(:, 3), '-b');
% xlabel('Sample'); ylabel('Magnitude'); legend('PC1', 'PC2', 'PC3');
% axis tight

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the square-rooted eigenvalue spectrum
% subplot(2, 2, 3); plot([1:1:31], flipdim(sqrt(diag(eigval)),1), 'o', 'MarkerSize', 2);
% xlabel('Component number'); ylabel('sqrt(\lambda)');
% axis tight

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part (d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a scatter plot of PC1 score versus the PC2 score
% eigval = -eigval;
% eigvec = -eigvec;
%  
% xdif = bsxfun(@minus, Spikes, mean(Spikes,2));
% Z = xdif'*eigvec;    % Z contains the new coordinates (PC scores). In words, center high-dimensional data (subtract the mean), then project onto axis defined by eigenvectors  
% % xdif(:,1)'*v(:,1);   % to get each PC score 
% subplot(2, 2, 4); plot(score(:,1),score(:,2),'.','MarkerSize', 5)   % scatter plot in the 2-d space using the PC1 and PC2
% xlabel('PC1 score'); ylabel('PC2 score');
% return;

%% Problem 2
% pc12 = Z(:,end-1:end)';      % matrix consisting of PC1 and PC2
% [test, train] = kfolds(pc12, 4);      % create 4-fold cross-validation

pc12 = score(:,1:2)';
[test, train] = kfolds(pc12, 4);

plotCount = 0;
for k = 1:8   % the number of k (cluster)
    
    % initialize parameters
    tempInitParams.mu = InitParams.mu(:,1:k);   % initialize mu, depending on the number of class k
    
    for kk = 1:k    % the number of clusters
        tempInitParams.Sigma(:,:,kk) = InitParams.Sigma(:,:,1);     % initialize sigma, depending on the number of class k  
    end
        
    tempInitParams.pi = ones(1,k)./k;   % initialize pi 
    
    for f = 1:4      % the # of cross-validation (4-fold)
        %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (a) %%%%%%%%%%%%%%%%%%%%%%%%
        % Training
        [mu, sigma, ppi]=func_GMM(tempInitParams, train{f,1});    % train with EM algorithm
        % Test
        teD = test{f,1};    % test data
        telogMat = nan(k, size(test{1,1},2));       % cluster by data points
        const = -0.5 * 2 * log(2*pi);
        for kk = 1:k     % the # of classes (clusters)
            S = sigma(:,:,kk);
            xdif = bsxfun(@minus, teD, mu(:,kk));     % bsxfun applies element-by-element binary operation to two arrays with singleton expansion enabled
            term1 = -0.5 * sum((xdif' * inv(S)) .* xdif', 2); % N x 1
            term2 = const - 0.5 * log(det(S)) + log(ppi(kk)); % scalar
            telogMat(kk,:) = term1' + term2;   % log(N(Xn|MuK,SigmaK)) + log(pi)
        end
        
        % Evaluate log P({x}) for the test data
        astar = max(telogMat, [], 1);
        adif = bsxfun(@minus, telogMat, astar);
        tenLL = log(sum(exp(adif), 1)) + astar; % 1 x N   nLL = log likelihood, sum across the clusters
        teLL(f) = sum(tenLL);     % log likelihood, sum across the trials  
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (b) %%%%%%%%%%%%%%%%%%%%%%%%
        % Plot the Gaussian distribution from the first CV fold
        if (f==1)   % in case of the first CV fold
           plotCount = plotCount + 1;   % plot tracking
           figure(plotCount);   % generate a new figure
            % Plot the training data in blue
            plot(train{f,1}(1,:),train{f,1}(2,:),'b.');hold on;
            % Plot the training data in red
            plot(test{f,1}(1,:),test{f,1}(2,:),'r.'); hold on;
            for c=1:k   % # of clusters
                func_plotEllipse(mu(:,c),sigma(:,:,c));
                hold on;
            end
            xlabel('PC1 score'); ylabel('PC2 score');
            titleStr=sprintf('ClusterNum=%d',k);
            title(titleStr);    
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%% Part (c) %%%%%%%%%%%%%%%%%%%%%%%%
        % For K=3, plot the waveforms corresponding to GMM cluster centers
        if (k==3 && f==1)
           plotCount = plotCount + 1;   % plot tracking
           figure(plotCount);   % generate a new figure
           muwaveform = comp(:,1:2)*mu + repmat(mean(Spikes,2),1,3);    % same as 'comp(:,1:2)*mu'
           plot(muwaveform);    % plot the average waveforms
           % bsxfun(@plus, comp(:,1:2)*pc12, mean(Spikes,2));   % plot the entire waveforms by projecting the latent variable back to the original space           
        end
    end
    
    LLk(k) = sum(teLL);     % log likelihood, sum across the cross-validation folds

end

plotCount = plotCount + 1;   % plot tracking
figure(plotCount);   % generate a new figure
p = plot(LLk,'-*');




