function fig = scatter_row_by_row(collect_mat, ylimits)

fig = figure; 
hold on; 
for r = 1:size(collect_mat, 1) 
    x = rand(1, size(collect_mat, 2)).*0.1-0.05 + [1:size(collect_mat, 2)]; 
    plot(x, collect_mat(r, :), 'k:')
    scatter(x, collect_mat(r, :), 50, 'filled')
end
clearvars r


ylim(ylimits)

xlim([.5 3.5])
set(gca, 'TickDir', 'out')

avg_collect_mat = nanmean(collect_mat(sum(collect_mat>0, 2)==size(collect_mat, 2), :), 1); 
for rr = 1:size(collect_mat, 2)
    plot([rr-.2, rr+.2], avg_collect_mat(rr).*ones(1, 2), 'k', 'LineWidth', 2.5)
end

end