function hAx = rasterPlotCellSubplot(plotTitle, timeStampCells, subPlotSet, fullColorMap, alphaValue)
    %

    % Get the total number of cell arrays
    numCellArrays = length(timeStampCells);

    % Use 'cool' colormap and sample it
    %fullColormap = colormap('cool');
    cmapIndices = round(linspace(1, size(fullColorMap, 1), numCellArrays));
    cmap = fullColorMap(cmapIndices, :);
   
    % Create a subplot using the provided position
    hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
    hold on;

    currentTrialOffset = 0;

    % Loop through each cell array
    for c = 1:numCellArrays
        cellArray = timeStampCells{c};
        numTrials = length(cellArray);
        currentColor = cmap(c, :); % Get the color for the current cell array

        % Loop through each trial and plot the rasters
        for trial = 1:numTrials
            % Retrieve the time points for the current trial
            timePoints = cellArray{trial};

            % Plot a line for each time point in the current trial using the line function
            for tp = 1:length(timePoints)
                line([timePoints(tp), timePoints(tp)], [trial + currentTrialOffset - 0.5, trial + currentTrialOffset + 0.5], 'Color', [currentColor, alphaValue], 'LineWidth', 2);
            end
        end
        
        % Update the current trial offset
        currentTrialOffset = currentTrialOffset + numTrials;
    end

    % Set the y-axis to display trials from top to bottom
    title(plotTitle)
    set(gca, 'YDir', 'reverse', 'TickDir', 'out');
    xlim([floor(min(cell2mat(cellfun(@(x) min(cell2mat(x(:))), timeStampCells, 'UniformOutput', false)))), ceil(max(cell2mat(cellfun(@(x) max(cell2mat(x(:))), timeStampCells, 'UniformOutput', false))))]);
    ylim([0, currentTrialOffset + 1]);
    xlabel('Time (s)');
    ylabel('# Trial');
    axis tight
    hold off;
end
