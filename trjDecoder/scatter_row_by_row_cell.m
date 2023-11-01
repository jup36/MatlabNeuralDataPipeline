function scatter_row_by_row_cell(collect_cell, ylimits)

figure; hold on; 
for r = 1:size(collect_cell, 1) 
    x = rand(1, size(collect_cell, 2)).*0.1-0.05 + (1:size(collect_cell, 2)); 
    plot(x, cell2mat(collect_cell(r, :)), 'k:')
    scatter(x, cell2mat(collect_cell(r, :)), 50, 'filled')
end
clearvars r


ylim(ylimits)

xlim([.5 3.5])
set(gca, 'TickDir', 'out')

avg_collect_cell = nanmean(cell2mat(collect_cell)); 
for rr = 1:size(collect_cell, 2)
    plot([rr-.2, rr+.2], avg_collect_cell(rr).*ones(1, 2), 'k', 'LineWidth', 2.5)
end

end