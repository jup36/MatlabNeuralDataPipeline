%% Problem 7.1

clear all; close all; clear functions; clc

%%
addpath(genpath('/Volumes/RAID2/parkj/MATLAB'));

%% load signal and timestamps, modify these!
filedirectory = '/Volumes/RAID2/parkj/MATLAB/Neural_signal_processing_BY/PS7';
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
S = cov(Spikes');                   % covariance matrix (Input matrix needs to be observation-by-features)
[loadings, PCscore] = pca(Spikes'); % coeffs are eigenvalues; scores are eigenvectors; S*v = v*D   
% loadings: the weight by which each standardized original variable should be multiplied to get the component score

% Double checking the PC score 
centerWFs = Spikes-repmat(mean(Spikes,2),1,length(Spikes)); % centered waveforms (mean substracted)
wf1 = centerWFs(:,1);               % centered waveform1
wf1ScorePC1 = dot(wf1,loadings(:,1))/norm(loadings(:,1));


% Plot the first three principal components
% subplot(2, 2, 2); plot(eigvec(:, 1), '-r');
% hold on; plot(eigvec(:, 2), '-g'); hold on; plot(eigvec(:, 3), '-b');
% xlabel('Sample'); ylabel('Magnitude'); legend('PC1', 'PC2', 'PC3');
% axis tight

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the square-rooted eigenvalue spectrum
% subplot(2, 2, 3); plot([1:1:31], sqrt(eigval), 'o', 'MarkerSize', 2);
% xlabel('Component number'); ylabel('sqrt(\lambda)');
% axis tight

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part (d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a scatter plot of PC1 score versus the PC2 score
eigval = -eigval;
eigvec = -eigvec;
 
xdif = bsxfun(@minus, Spikes, mean(Spikes,2));
Z = xdif'*eigvec;    % Z contains the new coordinates (PC scores). In words, center high-dimensional data (subtract the mean), then project onto axis defined by eigenvectors  
% % xdif(:,1)'*v(:,1);   % to get each PC score 
% subplot(2, 2, 4); plot(Z(:,1),Z(:,2),'.','MarkerSize', 5)   % scatter plot in the 2-d space using the PC1 and PC2
% xlabel('PC1 score'); ylabel('PC2 score');
% return;

%% 
pc12 = Z(:,1:2)';      % matrix consisting of PC1 and PC2
[test, train] = kfolds(pc12,4);      % 4-fold cross-validation

for k = 1:8   % the number of k (cluster)
    
    % initialize parameters
    tempInitParams.mu = InitParams.mu(:,1:k);   % initialize mu, depending on the number of class k
    
    for kk = 1:k    % the number of clusters
        tempInitParams.Sigma(:,:,kk) = InitParams.Sigma(:,:,1);     % initialize sigma, depending on the number of class k  
    end
        
    tempInitParams.pi = ones(1,k)./k;   % initialize pi 
    
    [mu, Sigma, ppi, gam, LL] = func_GMM_crossval(tempInitParams, train, test);

end




