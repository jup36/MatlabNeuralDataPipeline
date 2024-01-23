function [histoRateMean, histoRateStd, histoRateSem] = rasterHistogramPlotCell(timeStampCells, binEdges)
    %
    binWidth = mode(diff(binEdges)); 

    % Get the total number of cell arrays
    numCellArrays = length(timeStampCells);

    histoC = cell(length(timeStampCells), 1); 

    % Loop through each cell array
    for c = 1:numCellArrays
        cellArray = timeStampCells{c};
        numTrials = length(cellArray);
        tempMat = []; 

        % Loop through each trial and plot the rasters
        for trial = 1:numTrials
            % Retrieve the time points for the current trial
            timePoints = cellArray{trial};

            hc = histcounts(timePoints, binEdges); 
            tempMat = [tempMat; hc]; 
        end

        histoC{c, 1} = tempMat;  
    end

    histoRateC = cellfun(@(a) smooth2a(a./binWidth, 0, 1), histoC, 'UniformOutput', false); 
    histoRateMean = cell2mat(cellfun(@mean, histoRateC, 'UniformOutput', false));  
    histoRateStd = cell2mat(cellfun(@std, histoRateC, 'UniformOutput', false)); 
    histoSqN = sqrt(cell2mat(cellfun(@(a) size(a, 1), histoC, 'UniformOutput', false))); 
    histoRateSem = histoRateStd./repmat(histoSqN, 1, size(histoRateStd, 2)); 

    plotMeanSem(histoRateMean, histoRateSem, -2:0.05:1.95); 

    % Set the y-axis to display trials from top to bottom
    set(gca, 'TickDir', 'out');
    xlim([floor(min(cell2mat(cellfun(@(x) min(cell2mat(x(:))), timeStampCells, 'UniformOutput', false)))), ceil(max(cell2mat(cellfun(@(x) max(cell2mat(x(:))), timeStampCells, 'UniformOutput', false))))]);
    xlabel('Time');
    ylabel('Trial #');
    title('Raster histogram Plot');
    hold off;
end
