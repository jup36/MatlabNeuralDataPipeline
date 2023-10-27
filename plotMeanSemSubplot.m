function hAx = plotMeanSemSubplot(plotTitle, meanData, semData, timeX, subPlotSet, fullColormap, groupNames)
    % Create a subplot using the provided position
    hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
    hold on;

    % Get the total number of cell arrays
    numbGroups = size(meanData, 1); 

    % Use 'cool' colormap and sample it
    cmapIndices = round(linspace(1, size(fullColormap, 1), numbGroups));
    cmap = fullColormap(cmapIndices, :);

    % Define the alpha level for the shaded SEM area
    alphaValue = 0.2;

    plotHandles = zeros(1, numbGroups); % Initialize an array to store plot handles

    for i = 1:numbGroups
        % Plot the mean data
        plotHandles(i) = plot(timeX, meanData(i, :), 'Color', cmap(i, :), 'LineWidth', 1.5);
        
        % Use fill to create the shaded SEM area
        fillX = [timeX, fliplr(timeX)];
        fillY = [meanData(i, :) + semData(i, :), fliplr(meanData(i, :) - semData(i, :))];
        fill(fillX, fillY, cmap(i, :), 'FaceAlpha', alphaValue, 'EdgeColor', 'none');
    end
    
    % Add legend
    if nargin == 7 % Check if groupNames is provided as an argument
        legend(plotHandles, groupNames, 'Location', 'best');
    end
    
    title(plotTitle)
    xlabel('Time (s)');
    ylabel('Z-score');
    set(gca, 'XTick', ceil(timeX(1)):1:floor(timeX(end)), 'TickDir', 'out');
    set(gca, 'XTickLabelRotation', 0);
    axis tight; % Set axis tight
    hold off;
end
