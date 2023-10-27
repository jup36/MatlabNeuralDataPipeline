function rasterPlotCell(cellArrays)
    % RASTERPLOTMULTI Generates a raster plot for multiple cell arrays
    % Each entry in the cellArrays is a cell array containing cells representing trials with timestamps.

    % Get the total number of cell arrays
    numCellArrays = length(cellArrays);

    % Generate a blue colormap
    cmap = [linspace(1, 0, numCellArrays)', linspace(1, 0, numCellArrays)', ones(numCellArrays, 1)];

    % Prepare a figure for the raster plot
    figure;
    hold on;

    currentTrialOffset = 0;

    % Loop through each cell array
    for c = 1:numCellArrays
        cellArray = cellArrays{c};
        numTrials = length(cellArray);
        currentColor = cmap(c, :); % Get the color for the current cell array

        % Loop through each trial and plot the rasters
        for trial = 1:numTrials
            % Retrieve the time points for the current trial
            timePoints = cellArray{trial};

            % Plot a line for each time point in the current trial using the line function
            for tp = 1:length(timePoints)
                line([timePoints(tp), timePoints(tp)], [trial + currentTrialOffset - 0.7, trial + currentTrialOffset + 0.7], 'Color', currentColor, 'LineWidth', 2);
            end
        end
        
        % Update the current trial offset
        currentTrialOffset = currentTrialOffset + numTrials;
    end

    % Set the y-axis to display trials from top to bottom
    set(gca, 'YDir', 'reverse');
    xlim([floor(min(cell2mat(cellfun(@(x) min(cell2mat(x)), cellArrays, 'UniformOutput', false)))), ceil(max(cell2mat(cellfun(@(x) max(cell2mat(x)), cellArrays, 'UniformOutput', false))))]);
    ylim([0, currentTrialOffset + 1]);
    xlabel('Time');
    ylabel('Trial #');
    title('Raster Plot');
    hold off;
end
