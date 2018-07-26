%Problem set 7 

clear all; close all; clear functions; clc

%%
% addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps, modify these!
filedirectory = 'C:\Users\jup36\Dropbox\NSP\PS7';    % directory for 'behavcell' data of all doses is saved 
cd(filedirectory)       
load('ps7_data');    % spike snippet data

S = cov(Spikes');    % covariance matrix (columns: random variables, rows: observations)
[V,D,~] = eig(S);       