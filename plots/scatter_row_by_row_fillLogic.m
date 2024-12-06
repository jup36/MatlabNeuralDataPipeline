function fig = scatter_row_by_row_fillLogic(collect_mat, ylimits, fill_logic)

fig = figure; hold on; 
for r = 1:size(collect_mat, 1) 
    x = rand(1, size(collect_mat, 2)).*0.1-0.05 + (1:size(collect_mat, 2)); 
    plot(x, collect_mat(r, :), 'k:')
    
    % Check if the current row should be filled based on fill_logic
    if fill_logic(r)
        s = scatter(x, collect_mat(r, :), 50, 'filled');
        s.AlphaData = ones(1, length(collect_mat(r, :)))./1.5;
        s.MarkerFaceAlpha = 'flat';
    else
        s = scatter(x, collect_mat(r, :), 50); 
    end
end
clearvars r

% Set limits and labels if provided
if exist('ylimits', 'var')
    ylim(ylimits)
end
%xlim([.5 2.5])
set(gca, 'TickDir', 'out')

% Calculate and plot the mean line
avg_collect_mat = nanmean(collect_mat, 1); 
for rr = 1:size(collect_mat, 2)
    plot([rr-.2, rr+.2], avg_collect_mat(rr).*ones(1, 2), 'k', 'LineWidth', 2.5)
end

end
