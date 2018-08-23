%This script is to perform VTA cell clustering based on their responses to
% the onset of the reward. The reward response is computed via ROC area under the curve
% , representing the neuronal discriminability between the baseline and the
% reward onset periods. PCA is run on the ROC data to extract the principal components 
% and PC scores. Using the scores of the top 3 principal components, EM
% algorithm is used for clustering. 

clear all; clear functions; clc

%% load signal and timestamps modify these! 
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));

filedirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT';
cd(filedirectory);
fileName = 'SUA_AGB_perievtFR_50msbin_VTA.mat';      % filename for the peri-evt FR
load(fileName,'dat','pellet','tstpbeh','base')   % load neuronal and behavioral data

saveDirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\VTA_cell_classification';
saveName = 'SUA_AGB_VTA_cellclass_50ms_step_092315';      

%% set user variables!
uv.step = 0.05;     % step in sec (50 ms)
uv.binsize = 0.1;    % binsize in sec (100 ms)  % default bin width used for elsewhere is 200 ms, 100 ms was used by Uchida
uv.edges = [-5:uv.step:5];     % 50 ms bins to be used for histc
uv.numbunit = size(dat,1);     % number of units
uv.time1 = -2;      % starting point of the peri-event time windows
uv.time2 = 2;       % ending point of the peri-event time windows

%% ROC analysis
cd('C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\PeriPellet_VTA_50step100bin_-500to1000')    % directory for the raster data, binned raster data
FileList = dir('*raster_data.mat');     % list of each raster data 50 ms steps moving from -500 to 1000 ms around the time of reward onset

ROCmat = zeros(uv.numbunit,1500/50);    % 1500 ms window binned in 50 ms

for u = 1:uv.numbunit   % # of units
    
    load(FileList(u).name,'binned_raster_data','raster_labels');    % load the current unit's binned data
    
    base100ms = sum(full(base.SC{u}(1:50,1:40)),2)/(2*10);  % 2 for 2 sec, 10 as 100 ms bin was used for the peri-pellet binned spike counts
    base100ms(:,2) = 0;
       
    for t = 1:size(binned_raster_data,2)    % # of timebins
         
        temp100ms = binned_raster_data(1:50,t); % Only include the block 1 trials' response for cell classification
        temp100ms(:,2) = 1;
                
        glmmat = [base100ms;temp100ms];
        glmmat(:,1) = glmmat(:,1)- mean(base100ms(:,1),1);
        
        [~,~,~,ROCmat(u,t)] = perfcurve(glmmat(:,2),glmmat(:,1),1); % get the ROC AUC values  
    end
    
    fprintf('completed unit #%d\n', u);     % report the progress (completion of each cross-validation sample)
end


%% pca 
[coeff,PCcomp] = pca(ROCmat); % pca; get eigen value and eigen vectors

xdif = bsxfun(@minus, ROCmat', mean(ROCmat,1)'); % center the data
Z = xdif'*PCcomp; % get PC scores

pc123 = Z(:,1:3); % take PC 1,2,3

%% Visualization
% peri-evt raster for the GLM variable of interest 
unitnumb = 3;       % the # of unit of interest (get an idea by looking coefficients)
tempSTcell = pellet.ST{unitnumb,1};       % STcell of the unit 
tempFRmat = pellet.FR{unitnumb,1};
 
% grouping variable
label = zeros(size(tempFRmat,1),1);
label(1:50)=1; label(51:100)=2; label(101:150)=3;   % block ID
 
perievtrasterblocked( tempSTcell, uv.time1, uv.time2, label, 2 );   % peri-event raster plot, labelorder must be either -2 or 2! (Default = 2, -2 will flip the trial order when plotting)
set(gca,'XTick',[-1:1:1],'YTick',[50:50:200],'FontSize',14);
set(gca,'TickDir','out')
xlabel('Time (sec)','FontSize',14)
ylabel('Trials','FontSize',14)  

% % peri-evt filled line plots mean +- sem 
% perievtFRboundedline(tempFRmat, uv, label)
% set(gca,'XTick',[-1:1:1],'YTick',[0:5:200],'FontSize',14);
% xlabel('Time (sec)','FontSize',14)
% ylabel('Firing rate (Hz)','FontSize',14)  
% ylim([0 50])

figure;
xROC = [-0.5:0.05:0.95];
plot(xROC,ROCmat(unitnumb,:))
set(gca,'XTick',[-0.5:0.5:1],'YTick',[0:0.2:1],'FontSize',14);
set(gca,'TickDir','out')
ylim([0 1])
xlabel('Time (sec)','FontSize',14)
ylabel('AUC','FontSize',14)  

figure;
hold on;
plot3(pc123(:,1),pc123(:,2),pc123(:,3),'o')  % scatter plot
grid on
plot3(pc123(unitnumb,1),pc123(unitnumb,2),pc123(unitnumb,3),'ro')  

%% EM algorithm
% Initialize parameters
InitParams.mu = zeros(3,2);      % Initialization of the cluster center features(3) x clusters(2) 
InitParams.Sigma = zeros(3,3);   % Initialization of the covariance matrix. Assume that all cluster covariances are initialized to the same covariance matrix
InitParams.pi = [1/2, 1/2];      % Initialization of the prior cluster probability 

% cluster1
InitParams.mu(1,1) = 1.2;  % feature 1: amp
InitParams.mu(2,1) = 0;    % feature 2: halfwidth
InitParams.mu(3,1) = 0;    % feature 3: avFR

InitParams.Sigma(:,:,1) = diag(diag(cov(pc123)));   % Assume that all cluster covariances are initialized to the same covariance matrix 

% cluster2
InitParams.mu(1,2) = -1.5;   % feature 1
InitParams.mu(2,2) = 0;      % feature 2
InitParams.mu(3,2) = 0;      % feature 3

InitParams.Sigma(:,:,2) = diag(diag(cov(pc123)));   % Assume that all cluster covariances are initialized to the same covariance matrix 

[EM.MU,EM.SIGMA,EM.PI,EM.GAMMA,EM.LL] = func_GMM(InitParams,pc123');

[EM.maxGAMMA,EM.c] = max(EM.GAMMA);

% scatter plot with clustering color-coded
hold on;
for u = 1:length(pc123)   % # of data points
    if EM.c(1,u) == 1 % class 1
    plot3(pc123(u,1),pc123(u,2),pc123(u,3),'o','MarkerSize',6,'MarkerFaceColor',[76 135 198]./255,'MarkerEdgeColor','none')  % scatter for class1 units
    %scatter3(pc123(u,1),pc123(u,2),pc123(u,3),'MarkerFaceColor',[76 135 198]./255,'MarkerEdgeColor',[76 135 198]./255)  % scatter for class1 units
    else % class 2 
    plot3(pc123(u,1),pc123(u,2),pc123(u,3),'o','MarkerSize',6,'MarkerFaceColor',[220 30 62]./255,'MarkerEdgeColor','none')  % scatter for class1 units
    %scatter3(pc123(u,1),pc123(u,2),pc123(u,3),'MarkerFaceColor',[220 30 62]./255,'MarkerEdgeColor',[220 30 62]./255)  % scatter for class2 units
    end
end
hold off;
grid on

% plot log-likelihood versus iteration number
plot(EM.LL)   
title('Log likelihood versus iteration number')
xlabel('iteration #')
ylabel('log likelihood')
plot(EM.LL(2:100))
xlabel('iteration #')
ylabel('log likelihood')

cd(saveDirectory)
save(saveName,'ROCmat','dat','pellet','EM','pc*')

%% Auxiliary figures 
DAROC = ROCmat(EM.c==1',:); % ROCmat for DA cells
nDAROC = ROCmat(EM.c==2',:); % ROCmat for non-DA cells

DAROCidx = nanmean(DAROC(:,11:20),2);   % mean center ROC to be used for unit sorting
nDAROCidx = nanmean(nDAROC(:,11:20),2); % mean center ROC to be used for unit sorting

srtDAROC = sortrows([DAROC, DAROCidx],size(DAROC,2)+1); % sort rows of the ROC mat 
srtDAROC = srtDAROC(:,1:end-1); % get rid of the index from the ROCmat

imagesc(flipdim(srtDAROC,1))
%colormap(jet)



