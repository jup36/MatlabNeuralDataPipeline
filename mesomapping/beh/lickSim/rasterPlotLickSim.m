function hfig = rasterPlotLickSim(lickSimCell, numTrs)
%lickSimCell contains lick timestamps sampled from the exponential
% distribution fitted to the data

numSamples = unique(cell2mat(cellfun(@length, lickSimCell, 'UniformOutput', false)));
trs = randperm(numSamples, numTrs);

lickC = cell(size(lickSimCell, 1), 1); % tertiles
for ts = 1:size(lickSimCell, 1)
    % choose block (early, middle, late)
    tsLicksC = cellfun(@(a) a(trs), lickSimCell(ts, :), 'UniformOutput', false);

    % concatenate epochs (pre-cue, cue, post-cue)
    tsLickTsC = cellfun(@(a, b, c) [a; b; c], tsLicksC{1}, tsLicksC{2}, tsLicksC{3}, 'UniformOutput', false); % concatenate epochs
    lickC{ts} = tsLickTsC(:);
end

hfig = rasterPlotCell(lickC); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function hfig = rasterPlotCell(cellArrays)
        % RASTERPLOTMULTI Generates a raster plot for multiple cell arrays
        % Each entry in the cellArrays is a cell array containing cells representing trials with timestamps.

        % Get the total number of cell arrays
        numCellArrays = length(cellArrays);

        % Use 'cool' colormap and sample it
        fullColormap = colormap('cool');
        cmapIndices = round(linspace(1, size(fullColormap, 1), numCellArrays));
        cmap = fullColormap(cmapIndices, :);

        % Prepare a figure for the raster plot
        hfig = figure;
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

end