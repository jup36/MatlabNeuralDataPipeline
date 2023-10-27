function hAx = sigDrezDprmSubplot(plotTitle, sigDrez, subPlotSet, fullColorMap)

hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
hold on;

% Create an empty array to hold the scatter plot handles
h = [];

legendEntries = {};

if isfield(sigDrez{1}, 'dprime')
    dprime = cell2mat(cellfun(@(a) a.dprime, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(dprime), dprime, 100, fullColorMap(6, :), 'filled'); % Append handle
    plot(1:length(dprime), dprime, 'Color', fullColorMap(6, :), 'LineStyle', ':');
    legendEntries{end+1} = 'dprime';
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
