function blockByBlock4ChunksMeanSemPlot(stat4C)
% stat4C = stat.sR4c; 
hold on; 
for b = 1:4 % just include blocks of the first cycle for now
    [mb(b,:),~,sb(b,:)] = meanstdsem(cell2mat(cellfun(@(a) a', stat4C(:,b), 'un', 0))); % block-by-block mean and sem 
    x = [1:4]+(b-1)*4+0.1; 
    shadedErrorBar(x,mb(b,:),sb(b,:))
    scatter(x,mb(b,:), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor','k', 'MarkerFaceAlpha',.7)
end
hold off; 
end