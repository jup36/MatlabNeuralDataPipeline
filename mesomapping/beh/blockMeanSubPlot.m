function hAx = blockMeanSubPlot(plotTitle, blockMeanC, subPlotSet, fullColorMap, groupNames)
% blockMeanC = {lat.rwdFstLatBlockMean; lat.pnsFstLatBlockMean};
% Use 'cool' colormap and sample it
%cmapIndices = round(linspace(1, size(fullColorMap, 1), length(blockMeanC)));
cmapIndices = round(linspace(1, size(fullColorMap, 1), length(groupNames)));
cmap = fullColorMap(cmapIndices, :);

hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
hold on;
scatterHandles = gobjects(length(blockMeanC), 1); % preallocate a handle array for scatter plots

for i = 1:length(blockMeanC)
    scatterHandles(i) = scatter(1:length(blockMeanC{i}), blockMeanC{i}, 100, cmap(i, :), 'filled');
    plot(1:length(blockMeanC{i}), blockMeanC{i}, 'Color', cmap(i, :), 'LineStyle', ':')
end

% Add the legend using the scatter plot handles and group names
legend(scatterHandles, groupNames, 'Location', 'best'); % 'best' will place the legend at the best location to avoid overlapping with the data.


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

title(plotTitle)
% Update xlim and ylim with the new ranges
xlim([minX - xMargin, maxX + xMargin]);
ylim([minY - yMargin, maxY + yMargin]);

% Set TickDir to 'out', YTick to have six ticks, and YTickLabels with one decimal point
xTicks = 1:length(blockMeanC{i});
yTicks = linspace(minY - yMargin, maxY + yMargin, 6);
yTickLabels = arrayfun(@(x) sprintf('%.2f', x), yTicks, 'UniformOutput', false);

set(gca, 'TickDir', 'out', 'XTick', xTicks, 'YTick', yTicks, 'YTickLabel', yTickLabels) 

hold off;
%print(h, fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101723', 'Figure', 'firstLickLatencyTertiles'), '-dpdf', '-vector');  % '-painters' ensures the output is vector graphics
end