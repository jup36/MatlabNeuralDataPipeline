%This script deals with the real neural data problem in the NSP problem set 3.

clear all; clear functions; clc

addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));

cd('C:\Users\jup36\Dropbox\NSP\PS3')

load('ps3_realdata')

%% training phase 
% Multivariate gaussian with shared covariance 
numbunit = 97;

% get the spike counts of the time periods of interest 
train_trial_sc = nan(size(train_trial,2),numbunit,size(train_trial,1));       % train matrix; reach (8) x neuron (97) x  x trial (91)
for i = 1:size(train_trial,1)     % number of trials (91)
    for j = 1:size(train_trial,2)     % number of class (8)
        train_trial_sc(j,:,i) = nansum(train_trial(i,j).spikes(:,351:550),2)';    % take the spike counts of the 200 ms period of interest (a single 200 ms bin starting 150 ms after the reach target turns on)    
    end
end
clearvars i j 

% MLE of parameters for the multivariate gaussian with shared covariance
prior = zeros(size(train_trial,2),1);   % reach (8) x 1
prior(:) = 1/8;   % MLE for p(Ck) prior 

mu = nanmean(train_trial_sc,3);     % reach (8) x neuron (97); trial-averaged

covmat = cell(size(train_trial,2),1);   % cell containing the covariance matrix for each class (reach(8) x 1)   
rs_train_trial_sc = cell(size(train_trial,2),1);
for i = 1:size(train_trial,2)   % the number of reach (8)
    rs_train_trial_sc{i,1} = reshape(train_trial_sc(i,:,:),[numbunit size(train_trial,1) 1])';     % shape the train_trial_sc matrix such that each column represents a vector of each unit containing all 91 trials in each class 
    covmat{i,1} = cov(rs_train_trial_sc{i,1});   % get the covariance matrix (97 by 97 matrix)   
    sharedcovmat(:,:,i) = prior(i,1).*covmat{i,1};   % prior probability times covariance matrix of each class to be summed below
end
clearvars i j 
sharedcov = sum(sharedcovmat,3);    % get the shared covovariance matrix (97 x 97), to take the diagonal use; diag(diag(sum(sharedcovmat,3)))

%% test phase problem 1 (gaussian pdf, shared covariance) and 2 (gaussian pdf, each class covariance) 
test_trial_sc = nan(size(test_trial,2),numbunit,size(test_trial,1));       % test matrix; reach (8) x neuron (97) x  x trial (91)

for i = 1:size(test_trial,1)     % number of trials (91)
    for j = 1:size(test_trial,2)     % number of class (8)
        test_trial_sc(j,:,i) = nansum(test_trial(i,j).spikes(:,351:550),2)';    % take the spike counts of the 200 ms period of interest (a single 200 ms bin starting 150 ms after the reach target turns on)  
    end
end
clearvars i j 

conp1 = cell(size(test_trial,1), size(test_trial,2));    % conditional probability cell; trial (91) x reach (8)  
conp2 = cell(size(test_trial,1), size(test_trial,2));    % conditional probability cell; trial (91) x reach (8)  

decod1 = nan(size(test_trial,1), size(test_trial,2));    % matrix that will contain decoding result (1: correct, 0: incorrect)
decod2 = nan(size(test_trial,1), size(test_trial,2));    % matrix that will contain decoding result (1: correct, 0: incorrect)

% this loop takes the test phase data (trial-by-trial, time period of interest) and decode using the parameters and pdf obtained using the training phase data   
for i = 1:size(test_trial,1)     % number of trials (91)
    for j = 1:size(test_trial,2)     % number of class (8) 
        % get the log(likelihood) for each class
        for u = 1:size(test_trial,2)    % number of class (8) 
            %conp{i,j}(1,u) = log(mvnpdf(test_trial_sc(j,:,i), mu(u,:), sharedcov)) + log(prior(u,:));   % likelihood using mvnpdf
            conp1{i,j}(1,u) = mu(u,:)*inv(sharedcov)*test_trial_sc(j,:,i)' - 1/2*(mu(u,:)*inv(sharedcov)*mu(u,:)') + log(prior(u,:));   % likelihood with the shared covariance 
            %conp2{i,j}(1,u) = mu(u,:)*inv(covmat{u,1})*test_trial_sc(j,:,i)' - 1/2*(mu(u,:)*inv(covmat{u,1})*mu(u,:)') + log(prior(u,:));     % likelihood with each class covariance       
        end
        clearvars u
        % was decoding correct?
        [~,I1] = max(conp1{i,j});   % find the maximum likelihood, which indicats the most probable class 
        decod1(i,j) = I1==j;   % decoding result (correct: 1, incorrect: 0) 
%       [~,I2] = max(conp2{i,j});   % find the maximum likelihood, which indicats the most probable class 
%       decod2(i,j) = I2==j;   % decoding result (correct: 1, incorrect: 0) 
    end
end
clearvars i j u 

p1decodrate = sum(decod1(:))./(size(decod1,1)*size(decod1,2));

%% test phase problem 3 (poisson pdf) using the 'poisspdf'
conp3 = cell(size(test_trial,1), size(test_trial,2));    % conditional probability cell; trial (91) x reach (8)  
decod3 = nan(size(test_trial,1), size(test_trial,2));    % matrix that will contain decoding result (1: correct, 0: incorrect)

for i = 1:size(test_trial,1)    % number of trials (91)
    for j = 1:size(test_trial,2)    % number of class (8)
        % get the log(likelihood) for each class
        for u = 1:size(test_trial,2)    % number of class (8) 
            conp3{i,j}(1,u) = sum(log(poisspdf(test_trial_sc(j,:,i),mu(u,:)))) + log(prior(u,:));   % likelihood 
        end
        % was decoding correct?
        [~,I3] = max(conp3{i,j});   % find the maximum likelihood, which indicats the most probable class 
        decod3(i,j) = I3==j;   % decoding result (correct: 1, incorrect: 0) 
    end
end
clearvars i j u 

p3decodrate = sum(decod3(:))./(size(decod3,1)*size(decod3,2));

%% test phase problem 3 (poisson pdf) using the pdf equation
mu(mu == 0) = 0.01;     % set the minimum variance to 0.01
for i = 1:size(test_trial,1)    % number of trials (91)
    for j = 1:size(test_trial,2)    % number of class (8)
        % get the log(likelihood) for each class
        % xmat = reshape(test_trial_sc(j,:,:),[97 91])';
        for u = 1:size(test_trial,2)    % number of class (8) 
            conp3{i,j}(1,u) = test_trial_sc(j,:,i)*log(mu(u,:))' - sum(log(factorial(test_trial_sc(j,:,i)))) - sum(mu(u,:)) + log(prior(u,:));
        end
        [~,I3] = max(conp3{i,j});   % find the maximum likelihood, which indicats the most probable class 
        decod3(i,j) = I3==j;   % decoding result (correct: 1, incorrect: 0)     
    end
end
clearvars i j u 

p3decodrate = sum(decod3(:))./(size(decod3,1)*size(decod3,2));
 








