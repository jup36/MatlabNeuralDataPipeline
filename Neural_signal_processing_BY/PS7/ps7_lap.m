%Problem set 7 

clear all; close all; clear functions; clc

% addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

% load signal and timestamps, modify these!
filedirectory = 'C:\Users\Junchol\Dropbox\NSP\PS7';    % directory for 'behavcell' data of all doses is saved 
cd(filedirectory)       
load('ps7_data');    % spike snippet data

% get the covariance matrix
S = cov(Spikes');    % covariance matrix (columns: random variables, rows: observations)
[v,D,~] = eig(S);    % eigen decomposition (S*v = v*D) 
%plot(v(:,29:31))   % eigenvector waveforms
%plot(flipdim(sqrt(diag(D)),1),'.')      % square-rooted eigenvalue spectrum

% get the PC scores
xdiff = bsxfun(@minus, Spikes, mean(Spikes,2));     % center the high-dimensional data
z = xdiff'*v;    % PC score (n by d). In words, center high-dimensional data, then project onto axis defined by eigenvectors
%plot(z(:,31),z(:,30),'.')   % a scatter plot of PC1 vs. PC2

