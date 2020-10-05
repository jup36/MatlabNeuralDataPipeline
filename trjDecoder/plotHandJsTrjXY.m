function plotHandJsTrjXY(plotC,jsXYtrj,colorTheme,figSavePath,figSaveName,numbTrjToPlot)
%plotC = xyTrjP1Intn{f,c}; 
%jsXYtrj = xyJsTrjPullAlignIntn{f,c}; 

plotC = -plotC; 
jsXYtrj = -jsXYtrj;

% get trajectories by the deviation from the median trajectory
medTrj = nanmedian(plotC,3); 
distTrj = squeeze(sum(sum((plotC-repmat(medTrj,[1,1,size(plotC,3)])).^2)));
distTrj(:,2) = 1:length(distTrj); 
srtDistTrj = sortrows(distTrj,1); 
trjToPlot = srtDistTrj(1:numbTrjToPlot,2); 

% get trajectories by the proximity to the joystick target position
distTrj1 = squeeze(min(sum((plotC-repmat(jsXYtrj(:,1,:),[1,size(plotC,2),1])).^2),[],2)); 
distTrj1(:,2) = 1:length(distTrj1); 
srtDistTrj1 = sortrows(distTrj1,1); 
trjToPlot1 = srtDistTrj1(1:numbTrjToPlot,2); 

%plotCC = plotC(:,:,trjToPlot); 
plotCC = plotC(:,:,trjToPlot1); 
jsCC = jsXYtrj(:,:,trjToPlot1); 

figure; hold on;
%scatter(jsXYtrj(1), jsXYtrj(2), 100, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0 0 0]) % draw joystick 
for j = 1:size(plotCC,3) % trial
    [cTheme] = TNC_CreateRBColormapJP(size(plotCC,2),colorTheme);
    x = plotCC(1,:,j);
    y = plotCC(2,:,j);
    c = cTheme(1:end,:);
    scatter(x(1), y(1), 50, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(end,:),'MarkerFaceAlpha',.2) % draw starting point hTrj
    scatter(x(end), y(end), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c(end,:),'MarkerFaceAlpha',.7) % draw last point hTrj
    patch([x nan],[y nan],[1:length(y) nan],'FaceColor','none','EdgeColor','interp','lineWidth',2)
    colormap(c)
    
    jsX = jsCC(1,:,j); 
    jsY = jsCC(2,:,j);
    %patch([jsX nan],[jsY nan],[1:length(jsY) nan],'FaceColor','none','EdgeColor','k','lineWidth',5)
    plot([jsX(1) jsX(end)], [jsY(1) jsY(end)],'k','lineWidth',2)
    
    xlim([-12 5])
    ylim([-5 15])
    pbaspect([1 1 1])
    set(gca,'tickDir','out')
end
print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
hold off;
end