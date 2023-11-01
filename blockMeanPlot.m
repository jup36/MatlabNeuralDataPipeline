function h = blockMeanPlot(blockMeanC)
% blockMeanC = {lat.rwdFstLatBlockMean; lat.pnsFstLatBlockMean};
% Use 'cool' colormap and sample it
fullColormap = colormap('cool');
cmapIndices = round(linspace(1, size(fullColormap, 1), length(blockMeanC)));
cmap = fullColormap(cmapIndices, :);

h = figure; hold on;
for i = 1:length(blockMeanC)
    scatter(1:length(blockMeanC{i}), blockMeanC{i}, 100, cmap(i, :), 'filled')
    plot(1:length(blockMeanC{i}), blockMeanC{i}, 'Color', cmap(i, :), 'LineStyle', ':')
end
% Determine the full range of x-values and y-values
allLengths = cellfun(@length, blockMeanC);
allValues = cell2mat(blockMeanC);

minX = 1;
maxX = max(allLengths);
minY = min(allValues(:));
maxY = max(allValues(:));

% Compute 10% margins for x and y axes
xMargin = 0.1 * (maxX - minX);
yMargin = 0.1 * (maxY - minY);

% Update xlim and ylim with the new ranges
xlim([minX - xMargin, maxX + xMargin]);
ylim([minY - yMargin, maxY + yMargin]);

% Set TickDir to 'out', YTick to have six ticks, and YTickLabels with one decimal point
xTicks = 1:length(blockMeanC{i});
yTicks = linspace(minY - yMargin, maxY + yMargin, 6);
yTickLabels = arrayfun(@(x) sprintf('%.1f', x), yTicks, 'UniformOutput', false);

set(gca, 'TickDir', 'out', 'XTick', xTicks, 'YTick', yTicks, 'YTickLabel', yTickLabels, 'XLabel', 'Blocks', 'YLabel', 'Latency for first lick (s)');
hold off;
%print(h, fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101723', 'Figure', 'firstLickLatencyTertiles'), '-dpdf', '-vector');  % '-painters' ensures the output is vector graphics
end