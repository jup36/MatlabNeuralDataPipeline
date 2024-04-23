function plotColorListWithNumbers(colorList)
    % Determine the number of colors
    numColors = size(colorList, 1);
    
    % Create figure
    figure;
    hold on;
    axis off; % Turn off the axis for a cleaner look
    
    % Set the height of each patch based on the number of colors
    patchHeight = 1 / numColors;
    
    for i = 1:numColors
        % Calculate the bottom y-coordinate of the patch
        yBottom = (i-1) * patchHeight;
        
        % Create a rectangle [x, y, width, height]
        rectangle('Position', [0, yBottom, 1, patchHeight], 'FaceColor', colorList(i,:), 'EdgeColor', 'none');
        
        % Add text at the center of the patch
        text(0.5, yBottom + patchHeight / 2, num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
    hold off;
end
