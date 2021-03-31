function trjMovieWithRawVideo_noTrjOnlyVideo(hTrjDir, movSaveDir, movSaveName)
%'trjMovie' takes a data matrix containing actual and estimated
% trajectories from cortical and striatal neural population activity
% (pltm), and generates a MPEG movies of those trajectories across time. 
% plotm is data matrix to be plotted: variable-by-time
%https://www.mathworks.com/help/matlab/import_export/convert-between-image-sequences-and-video.html

filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'; 
%filePath = 'S:\Junchol_Data\JS2p0\WR40_082019\Matfiles'; 
cd(filePath)
%hTrjDir = 'S:\Junchol_Data\JS2p0\WR40_082019\082019_WR40';
%hTrjDir = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/082019_WR40'; 
%movSaveDir = 'S:\Junchol_Data\JS2p0\WR40_082019\Matfiles\Figure\KalmanFilter_decoding'; 
%movSaveDir ='/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles/Figure/KalmanFilter_decoding';
vfList = dir(hTrjDir); 

%% get s.dat from preprocessKFdecode*
kfDir = dir('rezKFdecodeHTrjCtxStrVel*'); 
load(fullfile(kfDir(1).folder,kfDir(1).name),'s');

% load the preprocessing data
%kfPreDir = dir('preprocessKFdecodeHTrjCtxStr_*'); 
%preS = load(fullfile(kfPreDir.folder,kfPreDir.name),'s');
%preS = preS.('s'); 

valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtx, 'un', 0)); % valid trials
stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, s.dat.laserIdx, 'un', 0)); % stim trials

behDir = dir('jsTime1k_KinematicsTrajectories*'); 
load(fullfile(behDir(1).folder,fullfile(behDir(1).name)),'jkvt'); 

currentFolder = pwd;
clearvars F

ax = 2; % axis of interest to plot (X:1,Y:2,Z:3)  

%% preprocess trjC
%trjC = cell(1,sum(sum(valTrI&~stmTrI))); 
stimeJkvtTr = [s.time.jkvtTr]; 
cnt = 0; 
for c = 1:size(s.dat.state,2)
    for r = 1:size(s.dat.state,1)
        if valTrI(r,c) && ~stmTrI(r,c)
            cnt = cnt + 1; % count valid no-stim trial
            tJkvt = s.dat.trialJkvt{r,c}; % trial in jkvt
            tmpI = stimeJkvtTr==tJkvt; 
            tmpTrjT = s.time(tmpI).t1R:20:s.time(tmpI).tE; 
            % get trajectories 
            tmpCtx = s.dat.stateCtx{r,c}(ax,:); 
            tmpStr = s.dat.stateStr{r,c}(ax,:); 
            tmpAct = s.dat.state{r,c}(ax,:);      
            tmpPull = double(s.dat.pullIdx{r,c}(1,:)); 
            trjS(cnt).stimeJkvtTr = find(stimeJkvtTr==tJkvt); 
            trjS(cnt).trjT = tmpTrjT(1:length(tmpCtx));              
            % video frames 
            tmpPath = jkvt(tJkvt).sVideoInfo.path; % just use the SIDE view video here          
            vName1 = strfind(tmpPath,'\cam1'); 
            vfListI = cell2mat(cellfun(@(a) strcmpi(a,tmpPath(vName1+1:end)), {vfList.name},'un',0)); 
            trjS(cnt).videoPath = fullfile(vfList(vfListI).folder, vfList(vfListI).name); 
            trjS(cnt).frmT = jkvt(tJkvt).vFrameTime; % frame time            
            trjS(cnt).frmI = trjS(cnt).trjT(1)<=trjS(cnt).frmT & trjS(cnt).frmT<=trjS(cnt).trjT(end); % frame index 
            trjS(cnt).frmT1rTe = trjS(cnt).frmT(trjS(cnt).frmI); % frame time within the T1R to Te range
            % interpolation and smoothing
            %trjS(cnt).trjm = smooth2a(intm([tmpCtx; tmpStr; tmpAct; tmpPull], sum(trjS(cnt).frmI)),0,10); % cortex, striatum, actual hand trajectory 
            trjS(cnt).trjm(1:3,:) = smooth2a(intm([tmpCtx; tmpStr; tmpAct], sum(trjS(cnt).frmI)),0,10); % cortex, striatum, actual hand trajectory 
            trjS(cnt).trjm(4,:) = intm(tmpPull, sum(trjS(cnt).frmI))>.5; % pull index 
        end
    end
end
clearvars c r %jkvt

%% select trajectories with high Str-Act correlations
corrStrAct = cell2mat(cellfun(@(a) corr(a(2,:)',a(3,:)'), {trjS(:).trjm}, 'un', 0)); 
corrStrAct(2,:) = 1:length(corrStrAct); 
srtCorrStrAct = sortrows(corrStrAct',-1);
trjSS = trjS(srtCorrStrAct(1:10,2));

%% plot
colorMap = [[100 149 237]./255; [50 205 50]./255; [50 50 50]./255]; % colorMap for cortex and striatum
trjCtxStrAct = smooth2a(cell2mat({trjSS(:).trjm}),0,5);
trjCtxStrAct = -trjCtxStrAct(1:3,:); % FLIP the sign to make the trajectory direction to be more intuitive  
trjCtxStrAct(4,:) = cell2mat(cellfun(@(a) a(4,:), {trjSS(:).trjm}, 'un', 0)); 

hold on; grid off; 
trjLength = cell2mat(cellfun(@length, {trjSS(:).trjm}, 'un', 0))'; 
pYmn = min(trjCtxStrAct(:)); % for patching pull moments Y
pYmx = max(trjCtxStrAct(:)); % for patching pull moments Y
for j = 1:size(trjSS,2)
    if j == 1
        tmpX = 5; 
        tmpXpull = 0; 
    else
        tmpX = sum(trjLength(1:j-1))+5; 
        tmpXpull = sum(trjLength(1:j-1)); 
    end
    pX1 = find(trjSS(j).trjm(4,:),1,'first'); % for patching pull moments X
    pX2 = find(trjSS(j).trjm(4,:),1,'last');  % for patching pull moments X
    patch(tmpXpull+[pX1 pX2 pX2 pX1],[pYmn pYmn pYmx pYmx], [0.9290 0.6940 0.1250], 'FaceAlpha', .25, 'EdgeColor', 'none')
    scatter(tmpX,trjCtxStrAct(3,tmpX),30,'MarkerFaceColor',colorMap(3,:),'MarkerEdgeColor','none')
    scatter(tmpX,trjCtxStrAct(1,tmpX),30,'MarkerFaceColor',colorMap(1,:),'MarkerEdgeColor','none')
    scatter(tmpX,trjCtxStrAct(2,tmpX),30,'MarkerFaceColor',colorMap(2,:),'MarkerEdgeColor','none')
end
clearvars j 
plot(trjCtxStrAct(3,:),'Color',colorMap(3,:))
plot(trjCtxStrAct(1,:),'Color',colorMap(1,:))
plot(trjCtxStrAct(2,:),'Color',colorMap(2,:))
hold off;
set(gca,'tickDir','out')
axis tight
figName1 = 'representativeVelocityTracesCtxStrActual'; 
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics',figName1),'-dpdf','-painters','-bestfit')

%% generate a combined video 
hold on; grid on;

%title(movSaveName,'Interpreter','none'); % better not to add 
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ylim',[min(min(cell2mat({trjSS(:).trjm}))) max(max(cell2mat({trjSS(:).trjm})))])

%set(gca,'units','pixels','position',[100 100 1920 1080]); 

trjCtx = animatedline('LineWidth',2,'Color',colorMap(1,:),'MaximumNumPoints',300);
trjStr = animatedline('LineWidth',2,'Color',colorMap(2,:),'MaximumNumPoints',300);
trjAct = animatedline('LineWidth',2,'Color',colorMap(3,:),'MaximumNumPoints',300);

valPtTot = 0; 
for j = 1:size(trjSS,2)
    pltm = trjSS(j).trjm; % current trial's hand trajectory (cortex, striatum, actual)
    vId = VideoReader(fullfile(trjSS(j).videoPath));  
    fr = 0;
    valPtTr = 0; 
    while hasFrame(vId) && fr<=size(trjSS(j).frmT,2)
        fr = fr+1;
        img1 = readFrame(vId); % hTrj raw video
        %imshow(img1)
        if trjSS(j).frmI(fr)
            valPtTot = valPtTot+1; % valid frame that corresponds to a point in the trajectory
            valPtTr = valPtTr+1;
                % actual trj
                addpoints(trjAct, valPtTot, pltm(3,valPtTr));
                headAct = scatter(valPtTot,pltm(3,valPtTr),75,'filled','MarkerFaceColor',colorMap(3,:),'MarkerEdgeColor',colorMap(3,:));
                % cortex trj
                addpoints(trjCtx, valPtTot, pltm(1,valPtTr));
                headCtx = scatter(valPtTot,pltm(1,valPtTr),75,'filled','MarkerFaceColor',colorMap(1,:),'MarkerEdgeColor',colorMap(1,:));
                % striatum trj
                addpoints(trjStr, valPtTot, pltm(2,valPtTr));
                headStr = scatter(valPtTot,pltm(2,valPtTr),75,'filled','MarkerFaceColor',colorMap(2,:),'MarkerEdgeColor',colorMap(2,:));
                drawnow
                
                %set(gcf,'PaperPositionMode','auto')
                tightfig(gcf) % set path for tightfig 
                tmpImg2Name = ['pltm' sprintf('%03d',valPtTot) '.jpg']; 
                print(fullfile(movSaveDir,tmpImg2Name),'-dpng','-r600')
                img2 = imread(fullfile(movSaveDir,tmpImg2Name));
                if valPtTot == 1
                    rsHeight = size(img1,1); 
                    rsWidth = round(size(img2,2)*(size(img1,1)/size(img2,1))); 
                end
                rsImg2 = imresize(img2,[rsHeight,rsWidth]);
                img12 = [img1, rsImg2]; 
                hImg1 = figure;
                imshow(img1); %, rsImg2]); 
                F(valPtTot) = getframe(hImg1);              
                pause(0.01);
                delete(headCtx); delete(headStr); delete(headAct); close
        end
    end
end

cd(movSaveDir)
video = VideoWriter(movSaveName,'MPEG-4');
video.Quality = 100;
video.FrameRate = 20;
open(video)
writeVideo(video,F)
close(video)
cd(currentFolder)
endctxRchWvel_DepthBin(:,3)),1,0)) % ctx Z
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
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     