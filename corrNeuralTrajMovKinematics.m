function [] = corrNeuralTrajMovKinematics(filePath, fileNameNeural, fileNameBeh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fileNameNeural = 'IT01Str_121317_gpfa_reachStart_5D_25msBin.mat';
fileNameBeh = 'BehVariables.mat';

behFile = dir(fullfile(filePath,fileNameBeh));
load(behFile.name,'pos1','vel1')

neuralTrajFile = dir(fullfile(filePath,fileNameNeural));
load(neuralTrajFile.name,'gpfaResult')

xorthStack = [];
for t = 1:length(gpfaResult.seqTrain) % increment trial
    xorthStack(:,:,t)= gpfaResult.seqTrain(t).xorth; 
    maxPosStack(t,1) = max(pos1(gpfaResult.seqTrain(t).trialId, 1:1000),[],2);
    maxPosStack(t,2) = t;
    maxAbsVelStack(t,1) = max(abs(vel1(gpfaResult.seqTrain(t).trialId, 1:1000)),[],2);
    maxAbsVelStack(t,2) = t;
end
clearvars t

sortMaxPosStack = sortrows(maxPosStack,-1);
sortMaxAbsVelStack = sortrows(maxAbsVelStack,-1);

%avgXorthTraj = mean(xorthStack,3); % the trial-averaged neural population trajectories 
maxPosXorthTraj = mean(xorthStack(:,:,sortMaxPosStack(1:30,2)),3);
maxAbsVelXorthTraj = mean(xorthStack(:,:,sortMaxAbsVelStack(1:30,2)),3);

neuralTrajDistPos = zeros(size(xorthStack,3),size(xorthStack,2),size(xorthStack,1)+1); % dimension: trial x timeBin x xOrthDim+1
neuralTrajDistVel = zeros(size(xorthStack,3),size(xorthStack,2),size(xorthStack,1)+1); % dimension: trial x timeBin x xOrthDim+1
for t = 1:length(gpfaResult.seqTrain) % increment trial
    for dim = 1:size(xorthStack,1)+1
        if dim < size(xorthStack,1)+1
            neuralTrajDistPos(t,:,dim) = diag(pdist2(maxPosXorthTraj(dim,:)', gpfaResult.seqTrain(t).xorth(dim,:)','euclidean'))';
            neuralTrajDistVel(t,:,dim) = diag(pdist2(maxAbsVelXorthTraj(dim,:)', gpfaResult.seqTrain(t).xorth(dim,:)','euclidean'))';
        elseif dim == size(xorthStack,1)+1
            neuralTrajDistPos(t,:,dim) = diag(pdist2(maxPosXorthTraj(:,:)', gpfaResult.seqTrain(t).xorth(:,:)','euclidean'))';
            neuralTrajDistVel(t,:,dim) = diag(pdist2(maxAbsVelXorthTraj(:,:)', gpfaResult.seqTrain(t).xorth(:,:)','euclidean'))';
        end
    end
end

%% correlation 
rPos = zeros(size(xorthStack,2), size(xorthStack,1)+1); % dimension: timeBin x xOrthDim+1 
pPos = zeros(size(xorthStack,2), size(xorthStack,1)+1); % dimension: timeBin x xOrthDim+1 
rVel = zeros(size(xorthStack,2), size(xorthStack,1)+1); % dimension: timeBin x xOrthDim+1 
pVel = zeros(size(xorthStack,2), size(xorthStack,1)+1); % dimension: timeBin x xOrthDim+1 


for dim = 1:size(xorthStack,1)+1 % increment dimensions
    [rPos(:,dim),pPos(:,dim)] = corr(neuralTrajDistPos(:,:,dim),maxPosStack(:,1)); 
    [rVel(:,dim),pVel(:,dim)] = corr(neuralTrajDistVel(:,:,dim),maxAbsVelStack(:,1)); 
end

cmap = TNC_CreateRBColormap(100,'hot'); % generate a colormap for imagesc psth
imagescJP(flip(abs(rPos')),cmap,[0 0.5]); pbaspect([1 1 1]); % plot the trial-averaged z-scored PSTHs



%% plot trajectories 
highAmpTrials = sortMaxPosStack([1:2,4:6],2);
lowAmpTrials  = sortMaxPosStack([end-8, end-4:end-1],2);

figure; 
for i = 1:length(lowAmpTrials)
    dim1TrajHighAmp = smooth(gpfaResult.seqTrain(highAmpTrials(i)).xorth(1,:),5); 
    dim2TrajHighAmp = smooth(gpfaResult.seqTrain(highAmpTrials(i)).xorth(2,:),5); 
    dim3TrajHighAmp = smooth(gpfaResult.seqTrain(highAmpTrials(i)).xorth(3,:),5); 
    
    dim1TrajLowAmp = smooth(gpfaResult.seqTrain(lowAmpTrials(i)).xorth(1,:),5);
    dim2TrajLowAmp = smooth(gpfaResult.seqTrain(lowAmpTrials(i)).xorth(2,:),5);
    dim3TrajLowAmp = smooth(gpfaResult.seqTrain(lowAmpTrials(i)).xorth(3,:),5);
    
    hold on;
    plot3(dim1TrajHighAmp, dim2TrajHighAmp, dim3TrajHighAmp,'color', 'r','LineWidth',2)
    plot3(dim1TrajHighAmp(1,1), dim2TrajHighAmp(1,1), dim3TrajHighAmp(1,1), 'o','MarkerSize',10, 'MarkerFaceColor',[255,182,193]./255, 'MarkerEdgeColor',[255,182,193]./255)
    plot3(dim1TrajHighAmp(25,1), dim2TrajHighAmp(25,1), dim3TrajHighAmp(25,1), 'o','MarkerSize',10, 'MarkerFaceColor',[255,0,0]./255, 'MarkerEdgeColor',[255,0,0]./255)

    
    plot3(dim1TrajLowAmp, dim2TrajLowAmp, dim3TrajLowAmp,'color', 'b','LineWidth',2)
    plot3(dim1TrajLowAmp(1,1), dim2TrajLowAmp(1,1), dim3TrajLowAmp(1,1), 'o','MarkerSize',10, 'MarkerFaceColor',[135,206,250]./255, 'MarkerEdgeColor',[135,206,250]./255)
    plot3(dim1TrajLowAmp(25,1), dim2TrajLowAmp(25,1), dim3TrajLowAmp(25,1), 'o','MarkerSize',10, 'MarkerFaceColor',[0,0,255]./255, 'MarkerEdgeColor',[0,0,255]./255)

end
grid on; hold off;




