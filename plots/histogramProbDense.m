function [pdfValues, pdfPoints] = histogramProbDense(datCells, colorCells, alpha)
    % Validate input lengths
    if length(datCells) ~= length(colorCells)
        error('Angle and color cell arrays must be of the same length.');
    end
    
    figure; % Create a new figure for the histograms and PDFs
    hold on; % Hold on to plot multiple histograms and PDFs in the same figure
    
    % Iterate through each condition
    for i = 1:length(datCells)
        % Extract current condition's dat and color
        dat = datCells{i};
        color = colorCells{i};
        
        % Plot histogram
        %h = histogram(dat, 'Normalization', 'probability', 'FaceColor', color, 'FaceAlpha', alpha);
        
        % Calculate and plot PDF
        [pdfValues{1, i}, pdfPoints{1, i}] = ksdensity(dat, 'Support', 'positive');
        plot(pdfPoints{1, i}, pdfValues{1, i}, 'LineWidth', 2, 'Color', color);
    end
    
    hold off; % Release hold on the figure
    %xlabel('Angle (degrees)');
    ylabel('Probability');
    %title('Probability Density and Histogram of Reach Angles');
    %legend('Condition 1', 'PDF 1', 'Condition 2', 'PDF 2', ...); % Update legend as needed
end


% function histogramProbDense(datCells, colorCells, alpha, numBinCells)
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
%         [pdfValues, pdfPoints] = ksdensity(dat, 'BandWidth', 3);
%         plot(pdfPoints, pdfValues, 'LineWidth', 2, 'Color', color);
%     end
%     
%     hold off; % Release hold on the figure
%     xlabel('Angle (degrees)');
%     ylabel('Probability');
%     title('Probability Density and Histogram of Reach Angles');
%     legend('Condition 1', 'PDF 1', 'Condition 2', 'PDF 2'); % Update legend as needed or dynamically based on input
% end

