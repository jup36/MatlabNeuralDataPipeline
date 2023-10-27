function hFig = plotMeanSem(meanData, semData, timeX)

% Create figure and axes
hFig = figure; 
hold on;

% Get the total number of cell arrays
% (Assuming your timestamps are stored in a cell array called timeStampCells)
numbGroups = size(meanData, 1); 

% Use 'cool' colormap and sample it
fullColormap = colormap('cool');
cmapIndices = round(linspace(1, size(fullColormap, 1), numbGroups));
cmap = fullColormap(cmapIndices, :);

% Define the alpha level for the shaded SEM area
alphaValue = 0.2;

for i = 1:numbGroups
    % Plot the mean data
    plot(timeX, meanData(i, :), 'Color', cmap(i, :), 'LineWidth', 1.5);
    
    % Use fill to create the shaded SEM area
    fillX = [timeX, fliplr(timeX)];
    fillY = [meanData(i, :) + semData(i, :), fliplr(meanData(i, :) - semData(i, :))];
    fill(fillX, fillY, cmap(i, :), 'FaceAlpha', alphaValue, 'EdgeColor', 'none');
end

xlabel('Time');
ylabel('Value');
set(gca, 'XTick', ceil(timeX(1)):1:floor(timeX(end)))
hold off;


end
