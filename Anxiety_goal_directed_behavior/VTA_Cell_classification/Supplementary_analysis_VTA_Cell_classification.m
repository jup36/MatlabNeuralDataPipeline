%This script is to generate auxiliary figures on the dopaminergic cell
% classification. 

clear all; clear functions; clc

%% load signal and timestamps modify these! 
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));
fileDirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\VTA_cell_classification';
load('SUA_AGB_VTA_cellclass_50ms_step_092315'); 
load('VTA_SUA_unitClassification_120115','DAIdx');

%% Auxiliary figures AUC imagesc
% All VTA cells
ROCidx = nanmean(ROCmat(:,11:20),2);   % mean center ROC to be used for unit sorting
srtROC = sortrows([ROCmat, ROCidx],size(ROCmat,2)+1); % sort rows of the ROC mat 
srtROC = srtROC(:,1:end-1); % get rid of the index from the ROCmat

% DA cells AUC
DAROC = ROCmat(DAIdx,:); % ROCmat for DA cells
DAROCidx = nanmean(DAROC(:,11:20),2);   % mean center ROC to be used for unit sorting
srtDAROC = sortrows([DAROC, DAROCidx],size(DAROC,2)+1); % sort rows of the ROC mat 
srtDAROC = srtDAROC(:,1:end-1); % get rid of the index from the ROCmat

% Plot DA cells sorted AUC
figure;
imagesc(flipdim(srtDAROC,1))
caxis([-0.1 1.1])
set(gca,'TickDir','out')

% DA cells PC scores
DApc123 = pc123(DAIdx,:); % PC scores
DApc123idx = zeros(sum(DAIdx),1);
DApc123idx(32,1) = 1;
srtDApc123 = sortrows([DApc123,DAROCidx],size(DApc123,2)+1);
srtDApc123 = srtDApc123(:,1:end-1);

figure;
imagesc(flipdim(srtDApc123(~logical(DApc123idx),:),1))
set(gca,'TickDir','out')
colormap gray
caxis([-1.5 1.5])

% non-DA cells
nDAROC = ROCmat(~DAIdx,:); % ROCmat for non-DA cells
nDAROCidx = nanmean(nDAROC(:,11:20),2); % mean center ROC to be used for unit sorting
srtnDAROC = sortrows([nDAROC, nDAROCidx],size(nDAROC,2)+1); % sort rows of the ROC mat 
srtnDAROC = srtnDAROC(:,1:end-1); % get rid of the index from the ROCmat

% Plot non-DA cells sorted AUC
figure;
imagesc(flipdim(srtnDAROC,1))
caxis([-0.1 1.1])
set(gca,'TickDir','out')

% non-DA cells PC scores
nDApc123 = pc123(~DAIdx,:); % PC scores
srtnDApc123 = sortrows([nDApc123,nDAROCidx],size(nDApc123,2)+1);
srtnDApc123 = srtnDApc123(:,1:end-1);

figure;
imagesc(flipdim(srtnDApc123,1))
set(gca,'TickDir','out')
colormap gray
caxis([-1.5 1.5])