function plotCorrTwoGroups(rhoMat)
figure; 
% plot the max reach Amp
randX = -.1 + (.1+.1)*rand(100,1);
hold on
for i = 1:size(rhoMat,1)
    plot([1 2]+randX(i), rhoMat(i,:),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',15)
    plot([1 2]+randX(i), rhoMat(i,:),':')
end
hold off
xlim([0.7 2.3])
set(gca,'tickDir','out')
set(gca,'Ytick',0:0.05:1)
%ylim(ylimit)
end