function plotIdvDataPointsTwoGroups(datMat)
% plot the max reach Amp
hold on;
randX = -.1 + (.1+.1)*rand(100,1);
for b = 1:size(datMat,2) % block
    xOffset = (b-1)*2; 
    tempDat = squeeze(datMat(:,b,:));
    for f = 1:size(tempDat,1) % each file
        plot(xOffset+[1 2]+randX(f), tempDat(f,:),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',5)
        plot(xOffset+[1 2]+randX(f), tempDat(f,:),':')
    end
end
%hold off
%xlim([0.7 2.3])
set(gca,'tickDir','out')
%set(gca,'Ytick',0:0.05:1)
%ylim(ylimit)
end