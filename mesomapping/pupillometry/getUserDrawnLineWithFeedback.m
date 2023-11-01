function [point1, point2] = getUserDrawnLineWithFeedback(image)
% Create a figure to display the image
hFigure = figure;

% Display the image
imshow(image);

% Create an axes for drawing the line
hAxes = gca;

% Prompt the user to draw a line
title('Draw a line over the metal platform for head fixation (click and drag)');

% Initialize variables to store line coordinates
x_coords = [];
y_coords = [];

% Initialize a flag to indicate when the drawing has started
drawingStarted = false;

% Initialize a flag to indicate when the line drawing is complete
drawingComplete = false;

% Initialize the starting point
startingPoint = [];

% Create a callback function for mouse motion
set(hFigure, 'WindowButtonMotionFcn', @(src, event) drawLine(src, event, hAxes));

% Create a callback function for mouse press
set(hFigure, 'WindowButtonDownFcn', @(src, event) startDrawing(src, event));

% Create a callback function for mouse release
set(hFigure, 'WindowButtonUpFcn', @(src, event) endDrawing(src, event));

% Wait for the user to finish drawing the line
uiwait(hFigure);

% Close the figure
close(hFigure);

% Extract the coordinates of the two endpoints
if ~isempty(x_coords) && ~isempty(y_coords)
    point1 = [x_coords(1), y_coords(1)];
    point2 = [x_coords(end), y_coords(end)];
else
    point1 = [];
    point2 = [];
end

% Nested function to update the line in real-time
    function drawLine(src, event, axesHandle)
        if ~drawingComplete && drawingStarted
            % Get the current mouse position
            currentPoint = get(axesHandle, 'CurrentPoint');
            x = currentPoint(1, 1);
            y = currentPoint(1, 2);

            % Update the displayed image with the straight line
            imshow(image, 'Parent', axesHandle);
            hold on;
            plot(axesHandle, [startingPoint(1), x], [startingPoint(2), y], 'r', 'LineWidth', 2);
            hold off;
        end
    end

% Nested function to start the line drawing
    function startDrawing(src, event)
        drawingStarted = true;

        % Record the starting point
        startingPoint = get(hAxes, 'CurrentPoint');
        startingPoint = startingPoint(1, 1:2); % Only take x and y, ignore z
    end

% Nested function to end the line drawing
    function endDrawing(src, event)
        drawingComplete = true;

        % Record the ending point
        endingPoint = get(hAxes, 'CurrentPoint');
        endingPoint = endingPoint(1, 1:2); % Only take x and y, ignore z

        x_coords = [startingPoint(1), endingPoint(1)];
        y_coords = [startingPoint(2), endingPoint(2)];

        uiresume(src);
    end
end


