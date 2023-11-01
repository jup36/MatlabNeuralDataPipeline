function hAx = sigDrezSubplot(plotTitle, sigDrez, subPlotSet, fullColorMap)

hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
hold on;

% Create an empty array to hold the scatter plot handles
h = [];

legendEntries = {};

if isfield(sigDrez{1}, 'hitRate')
    hitRate = cell2mat(cellfun(@(a) a.hitRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(hitRate), hitRate, 100, fullColorMap(1, :), 'filled'); % Append handle
    plot(1:length(hitRate), hitRate, 'Color', fullColorMap(1, :), 'LineStyle', ':');
    legendEntries{end+1} = 'Hit';
end

if isfield(sigDrez{1}, 'missRate')
    missRate = cell2mat(cellfun(@(a) a.missRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(missRate), missRate, 100, fullColorMap(2, :), 'filled'); % Append handle
    plot(1:length(missRate), missRate, 'Color', fullColorMap(2, :), 'LineStyle', ':');
    legendEntries{end+1} = 'Miss';
end

if isfield(sigDrez{1}, 'FaRate')
    FaRate = cell2mat(cellfun(@(a) a.FaRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(FaRate), FaRate, 100, fullColorMap(3, :), 'filled'); % Append handle
    plot(1:length(FaRate), FaRate, 'Color', fullColorMap(3, :), 'LineStyle', ':');
    legendEntries{end+1} = 'FA';
end

if isfield(sigDrez{1}, 'CrRate')
    CrRate = cell2mat(cellfun(@(a) a.CrRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(CrRate), CrRate, 100, fullColorMap(4, :), 'filled'); % Append handle
    plot(1:length(CrRate), CrRate, 'Color', fullColorMap(4, :), 'LineStyle', ':');
    legendEntries{end+1} = 'CR';
end

xlim([0.9 length(sigDrez)+0.1]);
ylim([0 1]);

set(gca, 'TickDir', 'out');
xlabel('Blocks');
ylabel('Rate');
title(plotTitle);
hold off;

% Add legend using the scatter plot handles
if ~isempty(legendEntries)
    legend(h, legendEntries, 'Location', 'best');
end

end
