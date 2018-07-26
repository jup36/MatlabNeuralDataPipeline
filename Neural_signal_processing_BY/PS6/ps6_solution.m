
clear all; close all; clear functions; clc

%%
addpath(genpath('/Volumes/RAID2/parkj/MATLAB'));

%% load signal and timestamps, modify these!
filedirectory = '/Volumes/RAID2/parkj/MATLAB/Neural_signal_processing_BY/PS6';    % directory for 'behavcell' data of all doses is saved 
cd(filedirectory)       
load('ps6_data');    % spike snippet data 
f_0 = 30000;
K = 3;      % 3 classes
[D,N] = size(Spikes);       % Dimension (time points), Number of samples 
t_spike = (0:D-1)/f_0;
InitParams = InitParams1;

for k=1:K
 InitParams.Sigma(:,:,k) = InitParams.Sigma(:,:,1);
end
[MU,SIGMA,PI,GAMMA,LL] = func_GMM(InitParams,Spikes);
[maxGAMMA,c] = max(GAMMA);

subplot(2,1,1)
plot(LL)
title('Log likelihood versus iteration number')
xlabel('iteration #')
ylabel('log likelihood')
subplot(2,1,2)
plot(LL(2:10))
xlabel('iteration #')
ylabel('log likelihood')

figure
for k=1:K
subplot(K,1,k)
hold on
plot(t_spike,Spikes(:,c==k),'k');
plot(t_spike,MU(:,k),'r-','linewidth',2)
plot(t_spike,MU(:,k)+sqrt(diag(SIGMA(:,:,k))),'r--','linewidth',1.5)
plot(t_spike,MU(:,k)-sqrt(diag(SIGMA(:,:,k))),'r--','linewidth',1.5)
ylim([min(min(Spikes)) max(max(Spikes))])
xlabel('time (seconds)')
ylabel('potential (mV)');
title(sprintf('Cluster %i voltage versus time',k));
end
saveas(gcf,'ps6_sol_fig2.pdf');