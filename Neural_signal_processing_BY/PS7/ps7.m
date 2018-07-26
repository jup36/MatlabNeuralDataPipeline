
%Problem set 7 
clear all; close all; clear functions; clc

%%
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps, modify these!
filedirectory = 'C:\Users\jup36\Dropbox\NSP\PS7';    % directory for 'behavcell' data of all doses is saved 
cd(filedirectory)       
load('ps7_data');    % spike snippet data

%% Problem 1
% Plot the N raw spike snippets in a "voltage vs. time" plot
% plot(Spikes)

% Plot the eigenvector waveforms corresponding to the three largest eigenvalues
S = cov(Spikes');    % covariance matrix (columns: random variables, rows: observations)
[v,D] = eig(S);      % v is the right eigen vector (Sv = vD), D is diagonal eigen values 
% plot(v(:,29:31))

% Plot the square-rooted eigenvalue spectrum
% plot(flipdim(sqrt(diag(D)),1),'o')

% Create a scatter plot of the PC1 score vs. the PC2 score
xdif = bsxfun(@minus, Spikes, mean(Spikes,2));
Z = xdif'*v;    % Z contains the new coordinates (PC scores). In words, center high-dimensional data (subtract the mean), then project onto axis defined by eigenvectors  
% xdif(:,1)'*v(:,31);   % to get each PC score 
% plot(Z(:,31),Z(:,30),'.')   % scatter plot in the 2-d space using the PC1 and PC2

%% Problem 2
pc2 = Z(:,30:31)';      % matrix consisting of PC1 and PC2

[test, train] = kfolds(pc2,4);      % 4-fold cross-validation

for k = 1:8   % the number of k (cluster)
    
    % initialize parameters
    tempInitParams.mu = InitParams.mu(:,1:k);   % initialize mu, depending on the number of class k
    
    for kk = 1:k
        tempInitParams.Sigma(:,:,kk) = InitParams.Sigma(:,:,1);     % initialize sigma, depending on the number of class k  
    end
        
    tempInitParams.pi = ones(1,k)./k;   % initialize pi 
    
    [mu, Sigma, ppi, gam, LL] = func_GMM_crossval(tempInitParams, train, test);

end
