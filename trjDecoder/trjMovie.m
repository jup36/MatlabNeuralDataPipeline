function trjMovie(pltm, movSaveDir, movSaveName)
%'trjMovie' takes a data matrix containing actual and estimated
% trajectories from cortical and striatal neural population activity
% (pltm), and generates a MPEG movies of those trajectories across time. 
% plotm is data matrix to be plotted: variable-by-time
currentFolder = pwd;
clearvars F
colorMap = [[100 149 237]./255; [50 205 50]./255; [50 50 50]./255]; % colorMap for cortex and striatum
hold on; grid on;
title(movSaveName,'Interpreter','none');
set(gca,'xtick',[])

trjCtx = animatedline('LineWidth',2,'Color',colorMap(1,:),'MaximumNumPoints',300);
trjStr = animatedline('LineWidth',2,'Color',colorMap(2,:),'MaximumNumPoints',300);
trjAct = animatedline('LineWidth',2,'Color',colorMap(3,:),'MaximumNumPoints',300);

for ii = 1:size(pltm,2)
    % cortex trj
    addpoints(trjCtx, ii, pltm(1,ii));
    headCtx = scatter(ii,pltm(1,ii),75,'filled','MarkerFaceColor',colorMap(1,:),'MarkerEdgeColor',colorMap(1,:));
    % striatum trj
    addpoints(trjStr, ii, pltm(2,ii));
    headStr = scatter(ii,pltm(2,ii),75,'filled','MarkerFaceColor',colorMap(2,:),'MarkerEdgeColor',colorMap(2,:));
    % actual trj
    addpoints(trjAct, ii, pltm(3,ii));
    headAct = scatter(ii,pltm(3,ii),75,'filled','MarkerFaceColor',colorMap(3,:),'MarkerEdgeColor',colorMap(3,:));
    drawnow
    F(ii) = getframe(gcf);
    pause(0.01);
    delete(headCtx); delete(headStr); delete(headAct);
end
cd(movSaveDir)
video = VideoWriter(movSaveName,'MPEG-4');
video.Quality = 100;
open(video)
writeVideo(video,F)
close(video)
cd(currentFolder)
end