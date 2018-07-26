%% Problem 8
clear all; close all; clear functions; clc

%%
%addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps, modify these!
%filedirectory = 'C:\Users\jup36\Dropbox\NSP\PS8';    % labPC directory  
filedirectory = 'C:\Users\Junchol\Dropbox\NSP\PS8';    % labtop directory  
cd(filedirectory)       
load('ps8_data.mat');

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,p] = size(Xsim);
% Peform PCA
mu = mean(Xsim);
Xsim_centered = Xsim-repmat(mu,n,1);    % centering Xsim 
[U,D] = eig(Xsim_centered'*Xsim_centered);    % eigen decomposition
[eigenvals, indx] = sort(diag(D),'descend');    % sort the eigen values in a descending order
% arrange the evectors according to the magnitude of the eigenvalues
U = U(:,indx);  % reorder the eigenvectors
% find the projection of the data onto PC1
z_hat_PCA = Xsim_centered*U(:,1);   % PC score
%Xhat_PCA = z_hat_PCA*U(:,1)' + repmat(mu,[length(Xsim) 1]);

% dim_reduce_plot(Xsim,z_hat_PCA,U(:,1))

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[estParamsPPCA, PPCALL] = probpca(Xsim', 1);    % fit the probabilistic PCA

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_PPCA = estParamsPPCA.L*estParamsPPCA.L' + diag(estParamsPPCA.Ph);  % estimated covariance of X
Ezx_PPCA = (estParamsPPCA.L'*inv(C_PPCA)*Xsim_centered')';   % expected Z|X, 8x1 
Xhat_PPCA = Ezx_PPCA*estParamsPPCA.L' + repmat(estParamsPPCA.d',[length(Xsim) 1]);   % projection back to the high-d space, 8x2

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Uhat,Dhat,Vhat] = svd(estParamsPPCA.L);    % singular value decomposition
dim_reduce_plot(Xsim,Ezx_PPCA,estParamsPPCA.L)

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (e) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[estParamsFA, FALL] = fastfa(Xsim', 1);    % fit the FA

C_FA = estParamsFA.L*estParamsFA.L' + diag(estParamsFA.Ph); 
Ezx_FA = (estParamsFA.L'*inv(C_FA)*Xsim_centered')';   % expected Z|X, 8x1 
Xhat_FA = Ezx_FA*estParamsFA.L' + repmat(estParamsFA.d',[length(Xsim) 1]);   % projection back to the high-d space, 8x2

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (f) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim_reduce_plot(Xsim,Ezx_FA,estParamsFA.L)






