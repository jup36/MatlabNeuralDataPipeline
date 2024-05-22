function hAx = rasterPlotLickWaterCellGroups(plotTitle, lickTsCells, waterTsCells, fullColorMap, xWin, cutoffWin, alphaValue, groupNames)
    %
    % Check if the provided group names match the number of cell arrays
    if length(groupNames) ~= length(lickTsCells)
        error('Number of group names must match the number of cell groups');
    end

    % Get the total number of cell arrays
    numCellArrays = length(lickTsCells);

    % Sample colors from the fullColorMap
    cmapIndices = round(linspace(1, size(fullColorMap, 1), numCellArrays));
    cmap = fullColorMap(cmapIndices, :);
    waterColor = [0, 87, 217]./255; 

    % Create a subplot using the provided position
    hAx = figure; hold on;

    currentTrialOffset = 0;
    legendHandles = zeros(numCellArrays, 1); % Array to hold the handles for the legend

    % Loop through each cell array
    for c = 1:numCellArrays
        lickCellArray = lickTsCells{c};
        lickCellArray = cellfun(@(a) a(a>=cutoffWin(1) & a<=cutoffWin(end)), lickCellArray, 'UniformOutput', false); 
        waterCellArray = waterTsCells{c}; 
        waterCellArray = cellfun(@(a) a(a>=cutoffWin(1) & a<=cutoffWin(end)), waterCellArray, 'UniformOutput', false); 
        
        numTrials = length(lickCellArray);
        currentColor = cmap(c, :); % Get the color for the current cell array

        % Loop through each trial and plot the rasters
        for trial = 1:numTrials
            % Retrieve the time points for the current trial
            lickTs = lickCellArray{trial};
            waterTs = waterCellArray{trial};
            % Plot a line for each time point in the current trial using the line function
            for tp = 1:length(lickTs)
                h = line([lickTs(tp), lickTs(tp)], [trial + currentTrialOffset - 0.5, trial + currentTrialOffset + 0.5], 'Color', [currentColor, alphaValue], 'LineWidth', 2);
                if trial == 1 && tp == 1 % Capture the handle only once for the legend
                    legendHandles(c) = h;
                end
            end
            
            for tt = 1:length(waterTs)
               scatter(waterTs(tt), trial + currentTrialOffset, 15, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', waterColor, 'MarkerFaceAlpha', 0.5); 
            end
        end
        
        % Update the current trial offset
        currentTrialOffset = currentTrialOffset + numTrials;
    end

    % Set the y-axis to display trials from top to bottom
    title(plotTitle)
    set(gca, 'YDir', 'reverse', 'TickDir', 'out');
    xlim(xWin)
    ylim([0, currentTrialOffset + 1]);
    xlabel('Time (s)');
    ylabel('# Trial');

    % Add the legend using provided group names
    legend(legendHandles, groupNames, 'Location', 'northwest', 'Interpreter', 'none');

    hold off;
end
