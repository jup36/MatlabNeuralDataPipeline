function hFig = firstLatencyPlotGngAcrossSession(plotTitle, latGoC, latNogoC, fullColorMap)

% sanity check
assert(size(latGoC, 1), size(latNogoC))

% Get the total number of cell arrays
numbSes = size(latGoC, 1);

hFig = figure;
hold on;

legendEntries = {};
scatterHandles = [];
xAnchor = 0;
for i = 1:size(latGoC, 1)
    latGo = latGoC{i}; 
    latNogo = latNogoC{i}; 

    % Go Create scatter plot and store the handle
    scatter(xAnchor+(1:length(latGo)), latGo, 100, fullColorMap(1, :), 'filled');
    plot(xAnchor+(1:length(latGo)), latGo, 'Color', fullColorMap(1, :), 'LineStyle', ':'); 
    
    % Nogo
    scatter(xAnchor+(1:length(latNogo)), latNogo, 100, fullColorMap(end, :), 'filled');
    plot(xAnchor+(1:length(latNogo)), latNogo, 'Color', fullColorMap(end, :), 'LineStyle', ':'); 

    xAnchor = xAnchor + length(latGo) + 1;
end

xlim([0 xAnchor]);
%ylim([0 max(3, max(latGo)+0.1)]);

set(gca, 'TickDir', 'out', 'XTick', 1:2:100);
xlabel('Blocks/Sessions');
ylabel('latency (s)');
title(plotTitle);
hold off;

end
