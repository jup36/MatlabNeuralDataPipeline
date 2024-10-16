function surfMatrix( filePath, matrixImg, matrixPval, cAxisVal, colorScheme, titleImg )
%This function takes a matrix, and generate a surf image of the matrix
% using the specified colormap 
alpha = 0.01; 

% get the colorMaps
[cMap] = TNC_CreateRBColormapJP(100,colorScheme);
[rbMap] = TNC_CreateRBColormapJP(100,'rb');

cd(fullfile(filePath,'Figure'))

hold on; plot(1:size(matrixImg,1), 1:size(matrixImg,2), ':r', 'LineWidth', 2)
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
plot(1:size(matrixImg,1), 1:size(matrixImg,2), ':r', 'LineWidth', 2)
%s=inputname(end); 
title(titleImg)
print(titleImg,'-dpdf')
close; 

figure; 
matrixSigImg = zeros(size(matrixImg)); % just get the portion of the rho matrix with significant p-values
matrixSigImg(matrixPval<alpha) = matrixImg(matrixPval<alpha); % just get the portion of the rho matrix with significant p-values 
hold on; plot(1:size(matrixSigImg,1), 1:size(matrixSigImg,2), ':k', 'LineWidth', 2)
imagesc(1:size(matrixSigImg,1),1:size(matrixSigImg,2),matrixSigImg);
colormap(rbMap)
colorbar
view(0,90)
axis tight
grid off
set(gca,'TickDir','out')
caxis(cAxisVal)
pbaspect([1 1 1]); % plot box aspect (ratio)
plot(1:size(matrixSigImg,1), 1:size(matrixSigImg,2), ':k', 'LineWidth', 2)
title(strcat(titleImg,'_SigPortion'),'Interpreter', 'none')
print(strcat(titleImg,'_SigPortion'),'-dpdf')
close; 

cd(fullfile(filePath))

end

