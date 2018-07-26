%% Problem 7.1
clear all; close all; clear functions; clc

%%
%addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps, modify these!
%filedirectory = 'C:\Users\jup36\Dropbox\NSP\PS7';    % labPC directory  
filedirectory = 'C:\Users\Junchol\Dropbox\NSP\PS8';    % labtop directory  
cd(filedirectory)       
load('ps8_data.mat');

%% %%%%%%%%%%%%%%%%%%%%%%% problem1 part (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply eigendecomposition
S = cov(Xplan);     % covariance matrix 
[eigvec, eigval] = eig(S);      % eigen desomposition
eigval = flipdim(diag(eigval),1);      % simply flip dim.
top3eigval = sum(eigval(1:3))/sum(eigval)*100;  % percent of variance explained by the first three principal components

[comp, score]=pca(Xplan);   % run PCA
%plot(sqrt(eigval),'-*')    % plot the square-rooted eigenvalue spectrum
%xlabel('PC1 score'); ylabel('PC2 score');

%% %%%%%%%%%%%%%%%%%%%%%%% problem1 part (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcs123 = score(:,1:3);  % first 3 pc scores
%plot3(pcs123(:,1),pcs123(:,2),pcs123(:,3),'.','MarkerSize',5)

[~,idx] = binfold(length(pcs123),8); % get the logical fold(reach angle) idx
cmap = hsv(8);  % define the colormap hsv
for i = 1:8     % # of folds
    plot3(pcs123(idx{i},1),pcs123(idx{i},2),pcs123(idx{i},3),'.','Color',cmap(i,:)); hold on;   
end
hold off;
xlabel('PC 1','fontsize',14)
ylabel('PC 2','fontsize',14)
zlabel('PC 3','fontsize',14)
title('Reach data projected into three-dimensional PC space','fontsize',14)

%% %%%%%%%%%%%%%%%%%%%%%%% problem1 part (c) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Um3 = comp(:,1:3);      % Um is D by M matrix, indicating the contribution of the dth neuron to the mth principal component
imagesc(Um3')
colormap jet
colorbar
xlabel('# of neurons','fontsize',14)
ylabel('PC components','fontsize',14)

%% %%%%%%%%%%%%%%%%%%%%%%% problem2 part (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,p] = size(Xsim);
% Peform PCA
mu = mean(Xsim);
Xsim_centered = Xsim-repmat(mu,n,1);    % centering Xsim 
[U,D] = eig(Xsim_centered'*Xsim_centered);  % eigen decomposition
[eigenvals, indx] = sort(diag(D),'descend');
% arrange the evectors according to the magnitude of the eigenvalues
U = U(:,indx);

figure;
plot(Xsim(:,1),Xsim(:,2),'.k'); hold on;    % plot each point    
plot(mean(Xsim(:,1),1),mean(Xsim(:,2),1),'.g');





