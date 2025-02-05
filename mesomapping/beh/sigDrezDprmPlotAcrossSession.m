function hFig = sigDrezDprmPlotAcrossSession(plotTitle, sigDrezC, fullColorMap)

% Get the total number of cell arrays
numbSes = size(sigDrezC, 1);

% Use 'cool' colormap and sample it
cmapIndices = round(linspace(5, size(fullColorMap, 1), numbSes));
cmap = fullColorMap(cmapIndices, :);

hFig = figure;
hold on;

legendEntries = {};
scatterHandles = [];
xAnchor = 0;
for i = 1:size(sigDrezC, 1)
    dPrm = sigDrezC{i}; 

    % Create scatter plot and store the handle
    scatterHandle = scatter(xAnchor+(1:length(dPrm)), dPrm, 100, cmap(i, :), 'filled');
    scatterHandles = [scatterHandles, scatterHandle]; % Append handle
    plot(xAnchor+(1:length(dPrm)), dPrm, 'Color', cmap(i, :), 'LineStyle', ':');
    legendEntries{end+1} = sprintf('session #%d', i);
    xAnchor = xAnchor + length(dPrm);
end

xlim([0 xAnchor]);
%ylim([0 max(3, max(dPrm)+0.1)]);

set(gca, 'TickDir', 'out', 'XTick', 1:2:100);
xlabel('Blocks/Sessions');
ylabel('dPrime');
title(plotTitle);
pbaspect([3 1 1]); 
hold off;

% Add legend using only the scatter plot handles
if ~isempty(legendEntries)
    legend(scatterHandles, legendEntries, 'Location', 'best');
end

end
