function rasterPlot(cellArray, colorArray, alphaValue)
    % RASTERPLOT Generates a raster plot for the given cell array
    % Each cell in the cellArray represents a trial and contains time points
    % relative to the trial onset.
    % colorArray is an RGB array for setting the color of the rasters.

    % Determine the number of trials
    numTrials = length(cellArray);

    % Check if colorArray length matches cellArray length
    if size(colorArray, 1) == 1
        colorArray = repmat(colorArray, [numTrials, 1]); 
    elseif size(colorArray, 1) ~= numTrials
        error('The length of the color array should match the number of trials.');
    end

    % Prepare a figure for the raster plot
    figure;
    hold on;

    % Loop through each trial and plot the rasters
    for trial = 1:numTrials
        % Retrieve the time points for the current trial
        timePoints = cellArray{trial};
        currentColor = [colorArray(trial, :) alphaValue]; % Get the color for the current trial and append alpha value

        % Plot a line for each time point in the current trial using the line function
        for tp = 1:length(timePoints)
            line([timePoints(tp), timePoints(tp)], [trial-0.7, trial+0.7], 'Color', currentColor, 'LineWidth', 2);
        end
    end

    % Set the y-axis to display trials from top to bottom
    set(gca, 'YDir', 'reverse');
    xlim([floor(min(cell2mat(cellfun(@min, cellArray, 'UniformOutput', false)'))), ceil(max(cell2mat(cellfun(@max, cellArray, 'UniformOutput', false)')))]);
    ylim([0, numTrials + 1]);
    xlabel('Time');
    ylabel('Trial #');
    title('Raster Plot');
    hold off;
end


