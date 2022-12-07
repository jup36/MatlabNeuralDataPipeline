function decodedTrjPlotReachPull_new(movSaveDir, movSaveName)
%'trjMovie' takes a data matrix containing actual and estimated
% trajectories from cortical and striatal neural population activity
% (pltm), and generates a MPEG movies of those trajectories across time. 
% plotm is data matrix to be plotted: variable-by-time
%https://www.mathworks.com/help/matlab/import_export/convert-between-image-sequences-and-video.html

filePath = '/Volumes/dudmanlab/junchol/js2p0/WR40_082019/Matfiles'; 
cd(filePath)

%% get s.dat from KF decoding results for reach and pull trajectories 
kfDirR = dir('rezKFdecodeHTrjCtxStrPosVel_reach_new_*'); 
sR = load(fullfile(kfDirR.folder,kfDirR.name),'s');
sR = sR.('s'); 
valTrI = cell2mat(cellfun(@(a) ~isempty(a), sR.dat.spkCtxR, 'un', 0)); % valid trials
stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, sR.dat.laserIdxR, 'un', 0)); % stim trials
sR.dat.stateR(valTrI) = deal(cellfun(@(a) a(1:3, :), sR.dat.stateR(valTrI),'un',0)); 

med = nanmedian(cell2mat(cellfun(@(a) a(:,1), sR.dat.stateR(valTrI)','un',0)),2);
sR.dat.stateR(valTrI) = deal(cellfun(@(a) a-repmat(med,1,size(a,2)),sR.dat.stateR(valTrI),'un',0));

load(fullfile(kfDirR.folder, kfDirR.name),'rez_reach');
stateR_ctx = rez_reach{1, 7}.rpr_est_ctx;
stateR_str = rez_reach{1, 7}.rpr_est_str;
stateR_cg = rez_reach{1, 7}.rpr_est_cg;

kfDirP = dir('rezKFdecodeHTrjCtxStrPosVel_pull_new_*'); 
sP = load(fullfile(kfDirP.folder,kfDirP.name),'s');
sP = sP.('s'); 

sP.dat.stateP(valTrI) = deal(cellfun(@(a) a(1:3, :), sP.dat.stateP(valTrI),'un',0)); 

load(fullfile(kfDirP.folder, kfDirP.name),'rez_pull');
stateP_ctx = rez_pull{1, 7}.rpr_est_ctx;
stateP_str = rez_pull{1, 7}.rpr_est_str;
stateP_cg = rez_pull{1, 7}.rpr_est_cg;

medP1 = nanmedian(cell2mat(cellfun(@(a) a(:,1), sP.dat.stateP(valTrI)','un',0)),2);
sP.dat.stateP(valTrI) = cellfun(@(a) a-repmat(med,1,size(a,2)),sP.dat.stateP(valTrI),'un',0);

% get decoded trajectories with correction for coordinates that occurred
% due to the median subtraction 
stateP_ctx = cell(size(stateR_ctx,1), size(stateR_ctx,2)); 
stateP_ctx(valTrI) = deal(cellfun(@(a) a+repmat(medP1,1,size(a,2))-repmat(med,1,size(a,2)), rez_pull{1, 7}.rpr_est_ctx(valTrI), 'un', 0));

stateP_str = cell(size(stateR_str,1), size(stateR_str,2)); 
stateP_str(valTrI) = deal(cellfun(@(a) a+repmat(medP1,1,size(a,2))-repmat(med,1,size(a,2)), rez_pull{1, 7}.rpr_est_str(valTrI), 'un', 0));

stateP_cg = cell(size(stateR_cg,1), size(stateR_cg,2)); 
stateP_cg(valTrI) = deal(cellfun(@(a) a+repmat(medP1,1,size(a,2))-repmat(med,1,size(a,2)), rez_pull{1, 7}.rpr_est_cg(valTrI), 'un', 0));

currentFolder = pwd;
clearvars F

ax = 2; % axis of interest to plot (X:1,Y:2,Z:3)  

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
            tmpTrjT = sR.time(tmpI).t1R:20:sR.time(tmpI).tE; 
            %% get trajectories 
            % reach portion
            tmpCtxR = stateR_ctx{r,c}(ax,:); 
            tmpStrR = stateR_str{r,c}(ax,:); 
            tmpCgR = stateR_cg{r,c}(ax,:); 
            tmpActR = sR.dat.stateR{r,c}(ax,:);  
            % pull portion 
            tmpCtxP = stateP_ctx{r,c}(ax,:); 
            tmpStrP = stateP_str{r,c}(ax,:); 
            tmpCgP = stateP_cg{r,c}(ax,:); 
            tmpActP = sP.dat.stateP{r,c}(ax,:);      
                  
            trjS(cnt).sTimeJkvtTr = find(sTimeJkvtTr==tJkvt); 
            
            % interpolation and smoothing
            trjS(cnt).trjR = [tmpCtxR; tmpStrR; tmpCgR; tmpActR; zeros(1,length(tmpCtxR))]; % cortex, striatum, actual hand trajectory REACH portion 
            trjS(cnt).trjP = [tmpCtxP; tmpStrP; tmpCgP; tmpActP; ones(1,length(tmpCtxP))]; % cortex, striatum, actual hand trajectory REACH portion 
            trjS(cnt).trjRP = [trjS(cnt).trjR, trjS(cnt).trjP]; % cortex, striatum, actual hand trajectory REACH portion            
        end
    end
end
clearvars c r %jkvt

%% select trajectories with high Str-Act correlations
corrStrAct = cell2mat(cellfun(@(a) corr(a(2,:)',a(4,:)'), {trjS(:).trjRP}, 'un', 0)); 
corrStrAct(2,:) = 1:length(corrStrAct); 
%srtCorrStrAct = sortrows(corrStrAct',-1);
trjSS = trjS(srtCorrStrAct(1:15,2));

%% plot
colorMap = [[100 149 237]./255; [50 205 50]./255; [251 176 64]./255; [50 50 50]./255]; % colorMap for cortex and striatum
trjCtxStrCgAct = smooth2a(cell2mat({trjSS(:).trjRP}),0,2);
trjCtxStrCgAct = -trjCtxStrCgAct(1:4,:); % FLIP the sign to make the trajectory direction to be more intuitive  
trjCtxStrCgAct(5,:) = cell2mat(cellfun(@(a) a(5,:), {trjSS(:).trjRP}, 'un', 0)); 

hold on; grid off; 
trjLength = cell2mat(cellfun(@length, {trjSS(:).trjRP}, 'un', 0))'; 
pYmn = min(trjCtxStrCgAct(:)); % for patching pull moments Y
pYmx = max(trjCtxStrCgAct(:)); % for patching pull moments Y
for j = 1:size(trjSS,2)
    if j == 1
        tmpX = 4; 
        tmpXpull = 0; 
            
        plot(trjCtxStrCgAct(4,tmpX:trjLength(1)),'Color',colorMap(4,:))
        plot(trjCtxStrCgAct(1,tmpX:trjLength(1)),'Color',colorMap(1,:))
        plot(trjCtxStrCgAct(2,tmpX:trjLength(1)),'Color',colorMap(2,:))
        plot(trjCtxStrCgAct(3,tmpX:trjLength(1)),'Color',colorMap(3,:))
        
    else
        tmpX = sum(trjLength(1:j-1))+4; 
        tmpXpull = sum(trjLength(1:j-1)); 
    end
    pX1 = find(trjSS(j).trjRP(5,:),1,'first'); % for patching pull moments X
    pX2 = find(trjSS(j).trjRP(5,:),1,'last');  % for patching pull moments X
    patch(tmpXpull+[pX1 pX2 pX2 pX1],[pYmn pYmn pYmx pYmx], [0.9290 0.6940 0.1250], 'FaceAlpha', .25, 'EdgeColor', 'none')
    scatter(tmpX,trjCtxStrCgAct(4,tmpX),30,'MarkerFaceColor',colorMap(4,:),'MarkerEdgeColor','none')
    scatter(tmpX,trjCtxStrCgAct(1,tmpX),30,'MarkerFaceColor',colorMap(1,:),'MarkerEdgeColor','none')
    scatter(tmpX,trjCtxStrCgAct(2,tmpX),30,'MarkerFaceColor',colorMap(2,:),'MarkerEdgeColor','none')
    scatter(tmpX,trjCtxStrCgAct(3,tmpX),30,'MarkerFaceColor',colorMap(3,:),'MarkerEdgeColor','none')

    xx = tmpX:sum(trjLength(1:j-1))+trjLength(j); 
    plot(xx, trjCtxStrCgAct(4, xx),'Color',colorMap(4,:))
    plot(xx, trjCtxStrCgAct(1, xx),'Color',colorMap(1,:))
    plot(xx, trjCtxStrCgAct(2, xx),'Color',colorMap(2,:))
    plot(xx, trjCtxStrCgAct(3, xx),'Color',colorMap(3,:))

end
clearvars j 
hold off;
set(gca,'tickDir','out')
axis tight
figName1 = 'representativeVelocityTracesCtxStrCgActual_Y'; 
print(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData/kfDecodeCorrKinematics',figName1),'-dpdf','-painters','-bestfit')

end