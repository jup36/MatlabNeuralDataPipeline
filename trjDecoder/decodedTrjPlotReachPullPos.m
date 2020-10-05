function decodedTrjPlotReachPullPos(filePath)
%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'; 
%filePath = 'S:\Junchol_Data\JS2p0\WR40_082019\Matfiles'; 
cd(filePath)
%hTrjDir = 'S:\Junchol_Data\JS2p0\WR40_082019\082019_WR40';
%hTrjDir = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/082019_WR40'; 
%movSaveDir = 'S:\Junchol_Data\JS2p0\WR40_082019\Matfiles\Figure\KalmanFilter_decoding'; 
%movSaveDir ='/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles/Figure/KalmanFilter_decoding';
%vfList = dir(hTrjDir); 

%% get s.dat from KF decoding results for reach and pull trajectories 
preDir = dir('preprocessKFdecodeHTrjCtxStr_reachPull_hTrjF_*'); 
preS = load(fullfile(preDir.folder,preDir.name),'s'); 
preS = preS.('s'); 

kfDirR = dir('rezKFdecodeHTrjCtxStrPos_reach_hTrjF_*'); 
sR = load(fullfile(kfDirR.folder,kfDirR.name),'s');
sR = sR.('s'); 

kfDirP = dir('rezKFdecodeHTrjCtxStrPos_pull_hTrjF_*'); 
sP = load(fullfile(kfDirP.folder,kfDirP.name),'s');
sP = sP.('s'); 

valTrI = cell2mat(cellfun(@(a) ~isempty(a), sR.dat.spkCtxR, 'un', 0)); % valid trials
stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, sR.dat.laserIdxR, 'un', 0)); % stim trials

behDir = dir('jsTime1k_KinematicsTrajectories*'); 
%load(fullfile(behDir(1).folder,fullfile(behDir(1).name)),'jkvt'); 

currentFolder = pwd;
clearvars F

ax = 1; % axis of interest to plot (X:1,Y:2,Z:3)  

%% revert global median subtraction 
medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), reshape(preS.dat.stateP(valTrI),1,1,[]),'un',0)),3);
medR1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), reshape(preS.dat.stateR(valTrI),1,1,[]),'un',0)),3); 
medD = medR1-medP1; 

% this is to revert the global mean subtraction applied separately to reach and pull trj 
for c = 1:size(sP.dat.stateP,2)
    for r = 1:size(sP.dat.stateP,1) 
        if ~isempty(sP.dat.stateP{r,c})
            sP.dat.statePmdCrt{r,c} = sP.dat.stateP{r,c}-repmat(medD(1:3,1),1,size(sP.dat.stateP{r,c},2));
            sP.dat.statePCtxMdCrt{r,c} = sP.dat.statePCtx{r,c}-repmat(medD(1:3,1),1,size(sP.dat.stateP{r,c},2)); 
            sP.dat.statePStrMdCrt{r,c} = sP.dat.statePStr{r,c}-repmat(medD(1:3,1),1,size(sP.dat.stateP{r,c},2)); 
        end
    end
end
clearvars r c 

%% preprocess trjC
%trjC = cell(1,sum(sum(valTrI&~stmTrI))); 
sTimeJkvtTr = [sR.time.jkvtTr]; 
cnt = 0; 
for c = 1:size(sR.dat.stateR,2)
    for r = 1:size(sR.dat.stateR,1)
        if valTrI(r,c) && ~stmTrI(r,c)
            cnt = cnt + 1; % count valid no-stim trial
            tJkvt = sR.dat.trialJkvt{r,c}; % trial in jkvt
            tmpI = sTimeJkvtTr==tJkvt; 
            %% get trajectories 
            % reach portion
            tmpCtxR = sR.dat.stateRCtx{r,c}(ax,:); 
            tmpStrR = sR.dat.stateRStr{r,c}(ax,:); 
            tmpActR = sR.dat.stateR{r,c}(ax,:);          
            % pull portion 
            tmpCtxP = sP.dat.statePCtxMdCrt{r,c}(ax,:); 
            tmpStrP = sP.dat.statePStrMdCrt{r,c}(ax,:); 
            tmpActP = sP.dat.statePmdCrt{r,c}(ax,:);      
                  
            trjS(cnt).sTimeJkvtTr = find(sTimeJkvtTr==tJkvt);                    
            % interpolation and smoothing
            trjS(cnt).trjR = [tmpCtxR; tmpStrR; tmpActR; zeros(1,length(tmpCtxR))]; % cortex, striatum, actual hand trajectory REACH portion 
            trjS(cnt).trjP = [tmpCtxP; tmpStrP; tmpActP; ones(1,length(tmpCtxP))]; % cortex, striatum, actual hand trajectory REACH portion 
            trjS(cnt).trjRP = [[tmpCtxR; tmpStrR; tmpActR; zeros(1,length(tmpCtxR))],[tmpCtxP; tmpStrP; tmpActP; ones(1,length(tmpCtxP))]]; % cortex, striatum, actual hand trajectory REACH portion 
            
        end
    end
end
clearvars c r %jkvt

%% select trajectories with high Str-Act correlations
corrStrAct = cell2mat(cellfun(@(a) corr(a(2,:)',a(3,:)'), {trjS(:).trjRP}, 'un', 0)); 
corrStrAct(2,:) = 1:length(corrStrAct); 
srtCorrStrAct = sortrows(corrStrAct',-1);
trjSS = trjS(srtCorrStrAct(1:15,2));

%% plot
colorMap = [[100 149 237]./255; [50 205 50]./255; [50 50 50]./255]; % colorMap for cortex and striatum
trjCtxStrAct = smooth2a(cell2mat({trjSS(:).trjRP}),0,2);
trjCtxStrAct = -trjCtxStrAct(1:3,:); % FLIP the sign to make the trajectory direction to be more intuitive  
trjCtxStrAct(4,:) = cell2mat(cellfun(@(a) a(4,:), {trjSS(:).trjRP}, 'un', 0)); 

hold on; grid off; 
trjLength = cell2mat(cellfun(@length, {trjSS(:).trjRP}, 'un', 0))'; 
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
    pX1 = find(trjSS(j).trjRP(4,:),1,'first'); % for patching pull moments X
    pX2 = find(trjSS(j).trjRP(4,:),1,'last');  % for patching pull moments X
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
figName1 = 'representativeXPositionTracesCtxStrActual'; 
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics',figName1),'-dpdf','-painters','-bestfit')

end