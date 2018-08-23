function [s]=surfMatrix( matrixImg, cAxisVal, colorScheme )
%This function takes a matrix, and generate a surf image of the matrix
% using the specified colormap 

% get the colorMap
[cMap] = TNC_CreateRBColormapJP(100,colorScheme);

hold on; plot(1:size(matrixImg,1), 1:size(matrixImg,2), ':r')
imagesc(1:size(matrixImg,1),1:size(matrixImg,2),matrixImg);
%surf(1:size(matrixImg,1),1:size(matrixImg,2),matrixImg,'EdgeColor','none','LineStyle','none');
colormap(cMap)
colorbar
view(0,90)
axis tight
grid off
set(gca,'TickDir','out')
caxis(cAxisVal)
pbaspect([1 1 1]); % plot box aspect (ratio)
plot(1:size(matrixImg,1), 1:size(matrixImg,2), ':r')

s=inputname(1); 
title(s)

end

