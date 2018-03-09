%This script performs (old) joystick trajectory analysis. 
% It first detects reaches in non-stim trials, then also detects
% pseudo-reaches in stim trials. Then, it performs PCA analysis on detected
% reaches. PCA is run on non-stim reaches, PC scores for pseudo-reaches are
% gained by projecting to PCs of non-stim reaches. 

clear all; clear functions; clc;

%% Raw voltage traces from bin and meta files
addpath(genpath('/Volumes/RAID2/parkj/MATLAB')) % use the matlab folder on the network drive

fileDirectory = '/Volumes/RAID2/parkj/NeuralData/IT03_Ldms_M1_013118/Matfiles';  % must use the desktop hard drive for faster processing 
cd(fileDirectory)

load('BehVariables', 'positionData') % append the position/velocity data variables

%% get reach properties: non-stim reaches
[ reachStart, reachStop, reach0, pos1, pos2, xpos1, ypos1, xpos2, ypos2 ] = getReachTimesJP( positionData ); % all reach traces, aligned to start (pos1), to stop (pos2)
% decimate the pos1 trajectories
dcmPos1 = zeros(size(pos1,1),size(pos1,2)/10); % decimate pos1 by a factor of 10
for t = 1:size(pos1,1)
    dcmPos1(t,:) = decimate(pos1(t,:),10);
end
clearvars t 
% pca on the decimted pos1 trajectories 
[pos1Coeff,pos1Score] = pca(dcmPos1); % plot(pos1Coeff(:,2))

%% get pseudo-reach trajectories of the laser-stim trials
[ stmReachStart, stmReachMW, stmPos1, stmXpos1, stmYpos1 ] = getPseudoReachTimesForStimTrials( positionData, ts.stmLaser );
% decimate the stmPos1 trajectories
dcmStmPos1 = zeros(size(stmPos1,1),size(stmPos1,2)/10); % decimate pos1 by a factor of 10
for t = 1:size(stmPos1,1)
    dcmStmPos1(t,:) = decimate(stmPos1(t,:),10);
end
clearvars t 

stmPos1Score = dcmStmPos1*pos1Coeff; % projection of stmPos1 data to the principal components of the no-stim pos1 trajectories

hold on; plot(pos1Score(:,1),pos1Score(:,2),'.'); plot(stmPos1Score(:,1),stmPos1Score(:,2),'r.'); hold off;

%hold on; plot3(pos1Score(:,1),pos1Score(:,2),pos1Score(:,3),'.'); plot3(stmPos1Score(:,1),stmPos1Score(:,2),stmPos1Score(:,3),'r.'); hold off;


