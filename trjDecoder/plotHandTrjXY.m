function plotHandTrjXY(plotC,jsXY,colorTheme,figSavePath,figSaveName,numbTrjToPlot)
% plotC = xyTrjP1Intn(1,:);
plotC = -plotC; 
jsXY = -jsXY;
medTrj = nanmedian(plotC,3); 
distTrj = squeeze(sum(sum((plotC-repmat(medTrj,[1,1,size(plotC,3)])).^2)));
distTrj(:,2) = 1:length(distTrj); 
srtDistTrj = sortrows(distTrj,1); 
trjToPlot = srtDistTrj(1:numbTrjToPlot,2); 
plotCC = plotC(:,:,trjToPlot); 

figure; hold on;
scatter(jsXY(1), jsXY(2), 100, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0 0 0]) % draw joystick 
for j = 1:size(plotCC,3) % trial
    [cTheme] = TNC_CreateRBColormapJP(size(plotCC,2),colorTheme);
    x = plotCC(1,:,j);
    y = plotCC(2,:,j);
    c = cTheme(1:end,:);
    scatter(x(1), y(1), 50, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(end,:),'MarkerFaceAlpha',.2) % draw starting point hTrj
    scatter(x(end), y(end), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c(end,:),'MarkerFaceAlpha',.7) % draw last point hTrj
    patch([x nan],[y nan],[1:length(y) nan],'FaceColor','none','EdgeColor','interp','lineWidth',2)
    colormap(c)
    xlim([-12 5])
    ylim([-5 15])
    pbaspect([1 1 1])
    set(gca,'tickDir','out')
end
print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
hold off;
end