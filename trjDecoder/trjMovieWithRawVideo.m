function trjMovieWithRawVideo(hTrjDir, movSaveDir, movSaveName)
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
load(fullfile(kfDir.folder,kfDir.name),'s');

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
                hImg12 = figure;
                imshow([img1, rsImg2]); 
                F(valPtTot) = getframe(hImg12);              
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
end