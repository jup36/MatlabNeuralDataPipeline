function meanSemErrorBarOverTime(meanMat, semMat, x, colorMat)
    % Validate the input dimensions
    if size(meanMat, 1) ~= size(semMat, 1) || size(meanMat, 2) ~= size(semMat, 2)
        error('Mean matrix and SEM matrix must be of the same size.');
    end
    if size(colorMat, 1) ~= 4 || size(colorMat, 2) ~= 3
        error('Color matrix must be 4x3 with RGB values.');
    end

    % Plot each mean ± SEM with shaded error bars
    hold on;
    for i = 1:size(meanMat, 1)
        % Get the mean and SEM values for the current row
        meanValues = meanMat(i, :);
        semValues = semMat(i, :);

        % Get the color for the current row
        color = colorMat(i, :);

        % Plot the mean values
        plot(x, meanValues, 'Color', color, 'LineWidth', 1.5);

        % Create shaded error bars
        fill([x fliplr(x)], [meanValues + semValues fliplr(meanValues - semValues)], ...
             color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    hold off;

    % Set axis labels
    xlabel('Time');
    ylabel('Value');

    % Set grid and box
    grid on;
    box on;

    % Add legend
    %legend('show');

    % Add title
    %title('Mean ± SEM Over Time');
end
