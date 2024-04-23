function plotReachAnglesOnFig(figObj, angleCells, colorCells, meanColorCells, lineWidth)
    % Ensure figObj is a valid figure handle
    if ~ishandle(figObj) || ~strcmp(get(figObj, 'Type'), 'figure')
        error('The first input must be a valid figure handle.');
    end
    
    % Activate the specified figure without creating a new figure or altering its properties
    set(0, 'CurrentFigure', figObj);
    
    % Validate input lengths
    if length(angleCells) ~= length(colorCells) || length(angleCells) ~= length(meanColorCells)
        error('All input cell arrays must be of the same length.');
    end
    
    % Validate if lineWidth is provided, else set to default
    if nargin < 5 || isempty(lineWidth)
        lineWidth = 1; % Default line width
    end

    % Loop through each set of angles
    for i = 1:length(angleCells)
        angles_degrees = angleCells{i}; % Get current set of angles
        color = colorCells{i}; % Get corresponding color for individual vectors
        meanColor = meanColorCells{i}; % Get corresponding color for the mean vector

        % Convert angles to radians and rotate 90 degrees counterclockwise
        angles_radians = angles_degrees * (pi / 180) + pi/2;

        % Calculate vector components for individual vectors
        x_components = cos(angles_radians);
        y_components = sin(angles_radians);

        % Plot vectors for current condition
        quiver(zeros(size(x_components)), zeros(size(y_components)), x_components, y_components, 0, 'Color', color, 'LineWidth', lineWidth);
        
        % Calculate mean angle and rotate
        mean_angle_radians = atan2(mean(sin(angles_radians)), mean(cos(angles_radians)));

        % Calculate vector component for mean vector
        mean_x_component = 1.5*cos(mean_angle_radians);
        mean_y_component = 1.5*sin(mean_angle_radians);
        
        % Plot mean vector with twice the line width and without arrowhead using 'plot'
        plot([0, mean_x_component], [0, mean_y_component], 'Color', meanColor, 'LineWidth', lineWidth * 20);
    end
end
