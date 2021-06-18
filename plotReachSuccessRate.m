function plotReachSuccessRate(control, laser)

plotm = [control laser]; 
figure;
randX = -.1 + (.1+.1)*rand(200,1);
hold on
for i = 1:length(control)
    plot([1 2]+randX(i), [plotm(i,1) plotm(i,2)],'o','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor','k','MarkerSize',10); 
    plot([1 2]+randX(i), [plotm(i,1) plotm(i,2)],':'); 
end
hold off
xlim([0.7 2.3])
ylim([0 1.1])
set(gca,'tickDir','out')
set(gca,'Ytick',0:0.2:1)
pbaspect([1 1 1])
end