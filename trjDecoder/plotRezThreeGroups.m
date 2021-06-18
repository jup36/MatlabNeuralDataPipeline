function plotRezThreeGroups(rezMat)
% plot the max reach Amp
figure;
randX = -.1 + (.1+.1)*rand(100,1);
hold on
for i = 1:size(rezMat,1)
    plot([1 2 3]+randX(i), rezMat(i,:),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',15)
    plot([1 2 3]+randX(i), rezMat(i,:),':')
end
hold off
xlim([0.7 3.3])
ylim([max(0,min(rezMat(:))-0.05),max(rezMat(:))+0.05])
set(gca,'tickDir','out')
set(gca,'Ytick',0:0.1:1)
end
