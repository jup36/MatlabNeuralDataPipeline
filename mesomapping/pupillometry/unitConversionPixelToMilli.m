
% Specify the path to your video file
videoFilePath = 'your_video_file.mp4'; % Change this to your video file path

% Create a VideoReader object to read the video
videoReader = VideoReader(videoFilePath);

% Read the first frame
firstFrame = readFrame(videoReader);

function [point1, point2] = getUserDrawnLineWithFeedback(image)
    % Create a figure to display the image
    hFigure = figure;
    
    % Display the image
    imshow(image);
    
    % Create an axes for drawing the line
    hAxes = gca;
    
    % Prompt the user to draw a line
    title('Draw a line (click and drag)');
    
    % Initialize variables to store line coordinates
    x_coords = [];
    y_coords = [];
    
    % Initialize a flag to indicate when the line drawing is complete
    drawingComplete = false;
    
    % Create a callback function for mouse motion
    set(hFigure, 'WindowButtonMotionFcn', @(src, event) drawLine(src, event, hAxes));
    
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
        if ~drawingComplete
            % Get the current mouse position
            currentPoint = get(axesHandle, 'CurrentPoint');
            x = currentPoint(1, 1);
            y = currentPoint(1, 2);
            
            % Store the coordinates
            x_coords = [x_coords, x];
            y_coords = [y_coords, y];
            
            % Update the displayed image with the line
            imshow(image, 'Parent', axesHandle);
            hold on;
            plot(axesHandle, x_coords, y_coords, 'r', 'LineWidth', 2);
            hold off;
        end
    end

    % Nested function to end the line drawing
    function endDrawing(src, event)
        drawingComplete = true;
        uiresume(src);
    end
end


