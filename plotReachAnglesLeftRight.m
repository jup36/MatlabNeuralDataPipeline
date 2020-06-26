function plotReachAnglesLeftRight(anglesLeft, anglesRight)
figure;
randX = -.1 + (.1+.1)*rand(200,1);
hold on
for i = 1:length(anglesLeft)
    plot(1+randX(i), anglesLeft(i),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',10)
end
for i = 1:length(anglesRight)
    plot(2+randX(i), anglesRight(i),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',10)
end
hold off
xlim([0.7 2.3])
%ylim(ylimit)
set(gca,'tickDir','out')
%set(gca,'Ytick',0:0.05:1)
end