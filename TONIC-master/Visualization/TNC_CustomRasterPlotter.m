function [] = TNC_CustomRasterPlotter(xList,yList,height,width,color,alphaLev,clearFlag,figNum)

numEvents = numel(xList);

if clearFlag == 1
    figure(figNum);
    clf;
else
    hold on;
end

for i=1:numEvents
    currX = xList(i);
    currY = yList(i);
    patch([0 width width 0] + currX, [0 0 height height] + currY, 'k', 'FaceColor', color , 'FaceAlpha', alphaLev, 'LineStyle', 'none');
end
