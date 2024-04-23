function histogramProbDensity(datCells, colorCells, alpha, numBinCells)
    % Validate input lengths
    if length(datCells) ~= length(colorCells) || length(datCells) ~= length(numBinCells)
        error('All input cell arrays must be of the same length.');
    end
    
    figure; % Create a new figure for the histograms and PDFs
    
    % Iterate through each condition
    for i = 1:length(datCells)
        % Extract current condition's data and color
        dat = datCells{i};
        color = colorCells{i};
        
        % First y-axis for the histogram
        hold on;
        histogram(dat, 'Normalization', 'probability', 'FaceColor', color, 'FaceAlpha', alpha, 'NumBins', numBinCells{i});
        ylabel('Probability');
    end

    xlabel('Variable'); % Set this to your variable name
    title('Histogram');
end


% function histogramProbDensity(datCells, colorCells, alpha, numBinCells)
%     % Validate input lengths
%     if length(datCells) ~= length(colorCells)
%         error('Angle and color cell arrays must be of the same length.');
%     end
%     
%     figure; % Create a new figure for the histograms and PDFs
%     hold on; % Hold on to plot multiple histograms and PDFs in the same figure
%     
%     % Iterate through each condition
%     for i = 1:length(datCells)
%         % Extract current condition's dat and color
%         dat = datCells{i};
%         color = colorCells{i};
%         
%         % Plot histogram with specified number of bins
%         %h = histogram(dat, 'Normalization', 'probability', 'FaceColor', color, 'FaceAlpha', alpha, 'NumBins', numBinCells{i});
%         
%         % Calculate and plot PDF
%         [pdfValues, pdfPoints] = ksdensity(dat, 'BandWidth', 40);
%         plot(pdfPoints, pdfValues, 'LineWidth', 2, 'Color', color);
%     end
%     
%     hold off; % Release hold on the figure
%     %xlabel('Angle (degrees)');
%     ylabel('Probability');
%     %title('Probability Density and Histogram of Reach Angles');
%     legend('Condition 1', 'PDF 1', 'Condition 2', 'PDF 2'); % Update legend as needed or dynamically based on input
% end
