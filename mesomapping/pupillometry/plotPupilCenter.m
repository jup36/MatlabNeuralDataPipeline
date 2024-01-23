function plotPupilCenter(pupilCenter, sizeRatio, rgbColor)
    if nargin < 2
        sizeRatio = 1.5; % Default size ratio
    end
    if nargin < 3
        rgbColor = [0, 0, 1]; % Default color (blue) if not provided
    end

    % Number of points
    numPoints = size(pupilCenter, 1);
    
    % Initialize figure
    figure;
    hold on;
    
    % Size and opacity parameters
    initialSize = 10; % Adjust as needed
    finalSize = initialSize * sizeRatio;
    initialAlpha = 0.01;
    finalAlpha = 0.85;

    % Loop through each point
    for i = 1:numPoints
        % Calculate current size and alpha
        currentSize = initialSize + (finalSize - initialSize) * (i - 1) / (numPoints - 1);
        currentAlpha = initialAlpha + (finalAlpha - initialAlpha) * (i - 1) / (numPoints - 1);

        % Plot circle with specified color and alpha
        scatter(pupilCenter(i, 1), pupilCenter(i, 2), currentSize, 'Filled', ...
                'MarkerFaceColor', rgbColor, 'MarkerFaceAlpha', currentAlpha, ...
                'MarkerEdgeAlpha', currentAlpha);
        
        % Draw line to next point with the same color and alpha (if not the last point)
        if i < numPoints
            nextAlpha = initialAlpha + (finalAlpha - initialAlpha) * i / (numPoints - 1);
            line([pupilCenter(i, 1), pupilCenter(i+1, 1)], [pupilCenter(i, 2), pupilCenter(i+1, 2)], ...
                 'Color', [rgbColor, nextAlpha], 'LineWidth', 1.5);
        end
    end
    
    hold off;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Pupil Center Movement Over Time');
    set(gca, 'TickDir', 'out')
    % print(fullfile(filePath, 'Figure', 'pupil_center_trial_1'), '-dpdf', '-vector')

end
