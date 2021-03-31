
filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};
figSavePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics';

%% load corrRez structure and collect data
% organize corr data
[ctxRchCorr, strRchCorr, ctxstrRchCorr, ctxPullCorr, strPullCorr, ctxstrPullCorr] = organizeCorrRezReachPull(filePath);
% organize r2 data
[ctxRchR2, strRchR2, ctxstrRchR2, ctxPullR2, strPullR2, ctxstrPullR2] = organizeR2RezReachPull(filePath);
% organize weights (C) data 
[ ctxRchW_posM, strRchW_posM, ctxRchW_velM, strRchW_velM, ctxPullW_posM, strPullW_posM, ctxPullW_velM, strPullW_velM ] = organizeKFdecoderWeights(filePath); 

%% save
%save(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics','trj_Corr_R2_kfWeight_Collect_reachPull_posVel')) % saved on 3/10/21
%load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics','trj_Corr_R2_kfWeight_Collect_reachPull_posVel'))

%% post-process weights data weights for X pos Ctx Str
% reach phase
[ctxRchWpos_DepthBin,~] = postprocessWeights(cell2mat(cellfun(@(a) a(:,1:end-1), ctxRchW_posM, 'un', 0))); % exclude cell Ids at the last column
[strRchWpos_DepthBin,~] = postprocessWeights(cell2mat(cellfun(@(a) a(:,1:end-1), strRchW_posM, 'un', 0))); % exclude cell Ids at the last column

figure; hold on; 
plot(ctxRchWpos_DepthBin(:,4),smooth2a(abs(ctxRchWpos_DepthBin(:,1)),1,0)) % ctx X
plot(strRchWpos_DepthBin(:,4),smooth2a(abs(strRchWpos_DepthBin(:,1)),1,0)) % str X    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_xPos_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_posM)); 
[strPullWpos_DepthBin,~] = postprocessWeights(cell2mat(strPullW_posM)); 

% weights for X pos pull phase Ctx Str
figure; hold on; 
plot(ctxPullWpos_DepthBin(:,4),smooth2a(abs(ctxPullWpos_DepthBin(:,1)),1,0)) % ctx X
plot(strPullWpos_DepthBin(:,4),smooth2a(abs(strPullWpos_DepthBin(:,1)),1,0)) % str X    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_xPos_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for Y pos Ctx Str
% reach phase
[ctxRchWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_posM)); 
[strRchWpos_DepthBin,~] = postprocessWeights(cell2mat(strRchW_posM)); 

figure; hold on; 
plot(ctxRchWpos_DepthBin(:,4),smooth2a(abs(ctxRchWpos_DepthBin(:,2)),1,0)) % ctx Y
plot(strRchWpos_DepthBin(:,4),smooth2a(abs(strRchWpos_DepthBin(:,2)),1,0)) % str Y    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_yPos_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_posM)); 
[strPullWpos_DepthBin,~] = postprocessWeights(cell2mat(strPullW_posM)); 

% weights for Y pos pull phase Ctx Str
figure; hold on; 
plot(ctxPullWpos_DepthBin(:,4),smooth2a(abs(ctxPullWpos_DepthBin(:,2)),1,0)) % ctx Y
plot(strPullWpos_DepthBin(:,4),smooth2a(abs(strPullWpos_DepthBin(:,2)),1,0)) % str Y    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_yPos_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for Z pos Ctx Str
% reach phase
[ctxRchWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_posM)); 
[strRchWpos_DepthBin,~] = postprocessWeights(cell2mat(strRchW_posM)); 

figure; hold on; 
plot(ctxRchWpos_DepthBin(:,4),smooth2a(abs(ctxRchWpos_DepthBin(:,3)),1,0)) % ctx Z
plot(strRchWpos_DepthBin(:,4),smooth2a(abs(strRchWpos_DepthBin(:,3)),1,0)) % str Z    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_zPos_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_posM)); 
[strPullWpos_DepthBin,~] = postprocessWeights(cell2mat(strPullW_posM)); 

% weights for Y pos pull phase Ctx Str
figure; hold on; 
plot(ctxPullWpos_DepthBin(:,4),smooth2a(abs(ctxPullWpos_DepthBin(:,3)),1,0)) % ctx Z
plot(strPullWpos_DepthBin(:,4),smooth2a(abs(strPullWpos_DepthBin(:,3)),1,0)) % str Z    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_zPos_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for X,Y,Z pos Ctx Str
% reach phase
[ctxRchWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_posM)); 
[strRchWpos_DepthBin,~] = postprocessWeights(cell2mat(strRchW_posM)); 

figure; hold on; 
plot(ctxRchWpos_DepthBin(:,4),smooth2a(nansum(abs(ctxRchWpos_DepthBin(:,1:3)),2),1,0)) % ctx XYZ
plot(strRchWpos_DepthBin(:,4),smooth2a(nansum(abs(strRchWpos_DepthBin(:,1:3)),2),1,0)) % str XYZ  
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_xyzPos_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWpos_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_posM)); 
[strPullWpos_DepthBin,~] = postprocessWeights(cell2mat(strPullW_posM)); 

% weights for Y pos pull phase Ctx Str
figure; hold on; 
plot(ctxPullWpos_DepthBin(:,4),smooth2a(nansum(abs(ctxPullWpos_DepthBin(:,1:3)),2),1,0)) % ctx XYZ
plot(strPullWpos_DepthBin(:,4),smooth2a(nansum(abs(strPullWpos_DepthBin(:,1:3)),2),1,0)) % str XYZ  
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_xyzPos_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for X vel Ctx Str
% reach phase
[ctxRchWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_velM)); 
[strRchWvel_DepthBin,~] = postprocessWeights(cell2mat(strRchW_velM)); 

figure; hold on; 
plot(ctxRchWvel_DepthBin(:,4),smooth2a(abs(ctxRchWvel_DepthBin(:,1)),1,0)) % ctx X
plot(strRchWvel_DepthBin(:,4),smooth2a(abs(strRchWvel_DepthBin(:,1)),1,0)) % str X    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_xVel_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_velM)); 
[strPullWvel_DepthBin,~] = postprocessWeights(cell2mat(strPullW_velM)); 

% weights for X vel pull phase Ctx Str
figure; hold on; 
plot(ctxPullWvel_DepthBin(:,4),smooth2a(abs(ctxPullWvel_DepthBin(:,1)),1,0)) % ctx X
plot(strPullWvel_DepthBin(:,4),smooth2a(abs(strPullWvel_DepthBin(:,1)),1,0)) % str X    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_xVel_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for Y vel Ctx Str
% reach phase
[ctxRchWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_velM)); 
[strRchWvel_DepthBin,~] = postprocessWeights(cell2mat(strRchW_velM)); 

figure; hold on; 
plot(ctxRchWvel_DepthBin(:,4),smooth2a(abs(ctxRchWvel_DepthBin(:,2)),1,0)) % ctx Y
plot(strRchWvel_DepthBin(:,4),smooth2a(abs(strRchWvel_DepthBin(:,2)),1,0)) % str Y    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_yVel_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_velM)); 
[strPullWvel_DepthBin,~] = postprocessWeights(cell2mat(strPullW_velM)); 

% weights for Y vel pull phase Ctx Str
figure; hold on; 
plot(ctxPullWvel_DepthBin(:,4),smooth2a(abs(ctxPullWvel_DepthBin(:,2)),1,0)) % ctx Y
plot(strPullWvel_DepthBin(:,4),smooth2a(abs(strPullWvel_DepthBin(:,2)),1,0)) % str Y    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_yVel_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for Z vel Ctx Str
% reach phase
[ctxRchWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_velM)); 
[strRchWvel_DepthBin,~] = postprocessWeights(cell2mat(strRchW_velM)); 

figure; hold on; 
plot(ctxRchWvel_DepthBin(:,4),smooth2a(abs(ctxRchWvel_DepthBin(:,3)),1,0)) % ctx Z
plot(strRchWvel_DepthBin(:,4),smooth2a(abs(strRchWvel_DepthBin(:,3)),1,0)) % str Z    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_zVel_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_velM)); 
[strPullWvel_DepthBin,~] = postprocessWeights(cell2mat(strPullW_velM)); 

% weights for Y vel pull phase Ctx Str
figure; hold on; 
plot(ctxPullWvel_DepthBin(:,4),smooth2a(abs(ctxPullWvel_DepthBin(:,3)),1,0)) % ctx Z
plot(strPullWvel_DepthBin(:,4),smooth2a(abs(strPullWvel_DepthBin(:,3)),1,0)) % str Z    
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_zVel_CtxStr'),'-dpdf','-bestfit','-painters')

%% post-process weights data weights for X,Y,Z vel Ctx Str
% reach phase
[ctxRchWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxRchW_velM)); 
[strRchWvel_DepthBin,~] = postprocessWeights(cell2mat(strRchW_velM)); 

figure; hold on; 
plot(ctxRchWvel_DepthBin(:,4),smooth2a(nansum(abs(ctxRchWvel_DepthBin(:,1:3)),2),1,0)) % ctx XYZ
plot(strRchWvel_DepthBin(:,4),smooth2a(nansum(abs(strRchWvel_DepthBin(:,1:3)),2),1,0)) % str XYZ  
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_reachPhase_xyzVel_CtxStr'),'-dpdf','-bestfit','-painters')

% pull phase
[ctxPullWvel_DepthBin,~] = postprocessWeights(cell2mat(ctxPullW_velM)); 
[strPullWvel_DepthBin,~] = postprocessWeights(cell2mat(strPullW_velM)); 

% weights for Y vel pull phase Ctx Str
figure; hold on; 
plot(ctxPullWvel_DepthBin(:,4),smooth2a(nansum(abs(ctxPullWvel_DepthBin(:,1:3)),2),1,0)) % ctx XYZ
plot(strPullWvel_DepthBin(:,4),smooth2a(nansum(abs(strPullWvel_DepthBin(:,1:3)),2),1,0)) % str XYZ  
hold off
set(gca,'TickDir','out')
xlim([0.1 4])
print(fullfile(figSavePath,'kfDecodeW_pullPhase_xyzVel_CtxStr'),'-dpdf','-bestfit','-painters')

%% plot MODIFY THIS PART!
% r2 Ctx vs Str vs CtxStr (stat: ctx vs. str) reach position
plotRezThreeGroups([cell2mat(ctxRchR2.posXyzToverall),cell2mat(strRchR2.posXyzToverall),cell2mat(ctxstrRchR2.posXyzToverall)]);
print(fullfile(figSavePath,'r2overall_reachPhase_posDecode'),'-dpdf','-bestfit','-painters')
[~,stats.r2rchPosCtxStrPval,~,stats.r2rchPosCtxStrStats] = ttest(cell2mat(ctxRchR2.posXyzToverall),cell2mat(strRchR2.posXyzToverall));

% r2 Ctx vs Str vs CtxStr (stat: ctx vs. str) pull position
plotRezThreeGroups([cell2mat(ctxPullR2.posXyzToverall),cell2mat(strPullR2.posXyzToverall),cell2mat(ctxstrPullR2.posXyzToverall)]);
print(fullfile(figSavePath,'r2overall_pullPhase_posDecode'),'-dpdf','-bestfit','-painters')
[~,stats.r2pullPosCtxStrPval,~,stats.r2pullPosCtxStrStats] = ttest(cell2mat(ctxPullR2.posXyzToverall),cell2mat(strPullR2.posXyzToverall));

% r2 Ctx vs Str vs CtxStr (stat: ctx vs. str) reach velocity
plotRezThreeGroups([cell2mat(ctxRchR2.velXyzToverall),cell2mat(strRchR2.velXyzToverall),cell2mat(ctxstrRchR2.velXyzToverall)]);
print(fullfile(figSavePath,'r2overall_reachPhase_velDecode'),'-dpdf','-bestfit','-painters')
[~,stats.r2rchVelCtxStrPval,~,stats.r2rchVelCtxStrStats] = ttest(cell2mat(ctxRchR2.velXyzToverall),cell2mat(strRchR2.velXyzToverall));

% r2 Ctx vs Str vs CtxStr (stat: ctx vs. str) pull velocity
plotRezThreeGroups([cell2mat(ctxPullR2.velXyzToverall),cell2mat(strPullR2.velXyzToverall),cell2mat(ctxstrPullR2.velXyzToverall)]);
print(fullfile(figSavePath,'r2overall_pullPhase_velDecode'),'-dpdf','-bestfit','-painters')
[~,stats.r2pullVelCtxStrPval,~,stats.r2pullVelCtxStrStats] = ttest(cell2mat(ctxPullR2.velXyzToverall),cell2mat(strPullR2.velXyzToverall));

%% save
save(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics','trj_Corr_R2_kfWeight_Collect_reachPull_posVel'),'stats','-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotRezThreeGroups(rezMat)
% plot the max reach Amp
figure;
randX = -.1 + (.1+.1)*rand(100,1);
hold on
for i = 1:size(rezMat,1)
    plot([1 2 3]+randX(i), rezMat(i,:),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',15)
    plot([1 2 3]+randX(i), rezMat(i,:),':')
end
hold off
xlim([0.7 3.3])
ylim([max(0,min(rezMat(:))-0.05),max(rezMat(:))+0.05])
set(gca,'tickDir','out')
set(gca,'Ytick',0:0.1:1)
end

function [wXyzDepthBinMean, wXyzDepthBinMed] = postprocessWeights(wXyzDepth) 
    % bin the x,y,z weights by depth bins
    depthInterval = 0.05; % 50 micron interval
    depthBins = min(wXyzDepth(:,4)):depthInterval:max(wXyzDepth(:,4))+.05; 
    
    [~,~,dI] = histcounts(wXyzDepth(:,4), depthBins); 
    
    wXyzDepthBinMean = nan(length(depthBins),5); 
    wXyzDepthBinMed = nan(length(depthBins),5); 
    
    for d = 1:length(wXyzDepthBinMean)
        wXyzDepthBinMean(d,:) = [nanmean(wXyzDepth(dI==d,:),1), sum(dI==d)]; 
        wXyzDepthBinMed(d,:)  = [nanmedian(wXyzDepth(dI==d,:),1), sum(dI==d)]; 
    end
end

function [ ctxRchW_posM, strRchW_posM, ctxRchW_velM, strRchW_velM, ctxPullW_posM, strPullW_posM, ctxPullW_velM, strPullW_velM ] = organizeKFdecoderWeights(filePath)

for f = 1:length(filePath)
    %% load files
    rDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_reach_*'));
    rS = load(fullfile(rDir.folder, rDir.name));
    pDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_pull_*'));
    pS = load(fullfile(pDir.folder, pDir.name));
    
    bscDir = dir(fullfile(filePath{f},'binSpkCountSTRCTXWR*'));
    stc = load(fullfile(bscDir(1).folder,bscDir(1).name),'spkTimesCell');
    depths = cell2mat(cellfun(@(a) a(2), stc.spkTimesCell(4,:),'un', 0));
    
    % reach 
    ctxRchW_posM{f,1} = collectAllWeights(rS.s.dat.params{1,7}.C_ctxVal,depths);
    strRchW_posM{f,1} = collectAllWeights(rS.s.dat.params{1,7}.C_strVal,depths);
    
    ctxRchW_velM{f,1} = collectAllWeights(rS.s.dat.params{1,8}.C_ctxVal,depths);
    strRchW_velM{f,1} = collectAllWeights(rS.s.dat.params{1,8}.C_strVal,depths);
    
    % pull
    ctxPullW_posM{f,1} = collectAllWeights(pS.s.dat.params{1,7}.C_ctxVal,depths);
    strPullW_posM{f,1} = collectAllWeights(pS.s.dat.params{1,7}.C_strVal,depths);
    
    ctxPullW_velM{f,1} = collectAllWeights(pS.s.dat.params{1,8}.C_ctxVal,depths);
    strPullW_velM{f,1} = collectAllWeights(pS.s.dat.params{1,8}.C_strVal,depths);
    
    % organize corr
    fprintf('processed file# %d\n', f)
end

    function [ meanW ] = collectAllWeights(weightMat,depths)
        nW = numel(weightMat);
        meanW = nanmean(cell2mat(reshape(weightMat,[1,1,nW])),3); % average across all trial resamples
        meanW(:,end+1) = depths'./1000; % add the depth in the last column to sort
        meanW = sortrows(meanW(~isnan(sum(meanW,2)),:),size(meanW,2)); % sort weights by depths
    end

end

function [ctxRchCorr, strRchCorr, ctxstrRchCorr, ctxPullCorr, strPullCorr, ctxstrPullCorr] = organizeCorrRezReachPull(filePath)

for f = 1:length(filePath) 
    %% load files
    rDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_reach_*')); 
    rS = load(fullfile(rDir.folder, rDir.name)); 
    pDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_pull_*')); 
    pS = load(fullfile(pDir.folder, pDir.name)); 
    % organize corr    
    for k = 1:length(rS.corrRez.ctx)
        switch k
            case 1
                % reach position X
                ctxRchCorr.posXyz{f,1}(1,1) = rS.corrRez.ctx{k}.all; 
                strRchCorr.posXyz{f,1}(1,1) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.posXyz{f,1}(1,1) = rS.corrRez.ctxstr{k}.all; 
                % pull position X
                ctxPullCorr.posXyz{f,1}(1,1) = pS.corrRez.ctx{k}.all; 
                strPullCorr.posXyz{f,1}(1,1) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.posXyz{f,1}(1,1) = pS.corrRez.ctxstr{k}.all; 
            case 2
                % reach position Y
                ctxRchCorr.posXyz{f,1}(1,2) = rS.corrRez.ctx{k}.all; 
                strRchCorr.posXyz{f,1}(1,2) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.posXyz{f,1}(1,2) = rS.corrRez.ctxstr{k}.all; 
                % pull position Y
                ctxPullCorr.posXyz{f,1}(1,2) = pS.corrRez.ctx{k}.all; 
                strPullCorr.posXyz{f,1}(1,2) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.posXyz{f,1}(1,2) = pS.corrRez.ctxstr{k}.all; 
            case 3
                % reach position Z
                ctxRchCorr.posXyz{f,1}(1,3) = rS.corrRez.ctx{k}.all; 
                strRchCorr.posXyz{f,1}(1,3) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.posXyz{f,1}(1,3) = rS.corrRez.ctxstr{k}.all; 
                % pull position Z
                ctxPullCorr.posXyz{f,1}(1,3) = pS.corrRez.ctx{k}.all; 
                strPullCorr.posXyz{f,1}(1,3) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.posXyz{f,1}(1,3) = pS.corrRez.ctxstr{k}.all; 
            case 4
                % reach velocity X
                ctxRchCorr.velXyz{f,1}(1,1) = rS.corrRez.ctx{k}.all; 
                strRchCorr.velXyz{f,1}(1,1) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.velXyz{f,1}(1,1) = rS.corrRez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullCorr.velXyz{f,1}(1,1) = pS.corrRez.ctx{k}.all; 
                strPullCorr.velXyz{f,1}(1,1) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.velXyz{f,1}(1,1) = pS.corrRez.ctxstr{k}.all; 
            case 5
                % reach velocity X
                ctxRchCorr.velXyz{f,1}(1,2) = rS.corrRez.ctx{k}.all; 
                strRchCorr.velXyz{f,1}(1,2) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.velXyz{f,1}(1,2) = rS.corrRez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullCorr.velXyz{f,1}(1,2) = pS.corrRez.ctx{k}.all; 
                strPullCorr.velXyz{f,1}(1,2) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.velXyz{f,1}(1,2) = pS.corrRez.ctxstr{k}.all; 
            case 6
                % reach velocity X
                ctxRchCorr.velXyz{f,1}(1,3) = rS.corrRez.ctx{k}.all; 
                strRchCorr.velXyz{f,1}(1,3) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.velXyz{f,1}(1,3) = rS.corrRez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullCorr.velXyz{f,1}(1,3) = pS.corrRez.ctx{k}.all; 
                strPullCorr.velXyz{f,1}(1,3) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.velXyz{f,1}(1,3) = pS.corrRez.ctxstr{k}.all; 
            case 7
                % reach position XYZ (fit together)
                ctxRchCorr.posXyzT{f,1} = diag(rS.corrRez.ctx{k}.all)'; 
                strRchCorr.posXyzT{f,1} = diag(rS.corrRez.str{k}.all)'; 
                ctxstrRchCorr.posXyzT{f,1} = diag(rS.corrRez.ctxstr{k}.all)'; 
                % pull position XYZ (fit together)
                ctxPullCorr.posXyzT{f,1} = diag(pS.corrRez.ctx{k}.all)'; 
                strPullCorr.posXyzT{f,1} = diag(pS.corrRez.str{k}.all)'; 
                ctxstrPullCorr.posXyzT{f,1} = diag(pS.corrRez.ctxstr{k}.all)';        
            case 8
                % reach velocity XYZ (fit together)
                ctxRchCorr.velXyzT{f,1} = diag(rS.corrRez.ctx{k}.all)'; 
                strRchCorr.velXyzT{f,1} = diag(rS.corrRez.str{k}.all)';
                ctxstrRchCorr.velXyzT{f,1} = diag(rS.corrRez.ctxstr{k}.all)'; 
                % pull velocity XYZ (fit together)
                ctxPullCorr.velXyzT{f,1} = diag(pS.corrRez.ctx{k}.all)'; 
                strPullCorr.velXyzT{f,1} = diag(pS.corrRez.str{k}.all)'; 
                ctxstrPullCorr.velXyzT{f,1} = diag(pS.corrRez.ctxstr{k}.all)';  
        end
    end
    % organize corr    
    fprintf('processed file# %d\n', f)
end    
  
