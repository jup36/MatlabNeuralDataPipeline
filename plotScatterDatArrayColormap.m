function plotScatterDatArrayColormap(datArray, fullcmap)
% varargin can provide names to be used for the figure legend

numDatPoints = size(datArray, 1); 
x = 1:numDatPoints; 

% Create figure and axes
%hFig = figure; 
%hold on;

cmapIndices = round(linspace(1, size(fullcmap, 1), numDatPoints));
cmap = fullcmap(cmapIndices, :);

% Array to store handles for the mean data plots
plot(x, datArray, 'Color', [0.4 0.4 0.4], 'LineStyle', ':')

for i = 1:numDatPoints
    % Plot the mean data
    datArrayHandles(i) = scatter(x(i), datArray(i, :), 150, cmap(i, :), 'filled');
end

xlim([min(x)-0.1 max(x)+0.1])
set(gca, 'TickDir', 'out')
%hold off;

end
