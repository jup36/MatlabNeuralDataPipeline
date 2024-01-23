function vertline(x, ylimit, varargin)
    % vertline: Draws vertical lines on a plot at specified x-coordinates.
    % Usage: vertline(x, ylimit, 'lineColor', [r g b], 'lineStyle', style, 'lineWidth', width, 'alpha', alphaValue)
    %   x - A vector of x-coordinates where the lines should be drawn
    %   ylimit - Two-element vector specifying y-limits for the lines
    %   varargin - Name-value pairs for line properties

    % Create an input parser
    p = inputParser;
    addParameter(p, 'lineColor', [0, 0, 0], @(x) (isvector(x) && length(x) == 3) || (isvector(x) && length(x) == 4));
    addParameter(p, 'lineStyle', '-', @ischar);
    addParameter(p, 'lineWidth', 1, @isnumeric);
    addParameter(p, 'alpha', 1, @(x) isnumeric(x) && x >= 0 && x <= 1); % Default alpha: 1 (opaque)

    % Parse the input arguments
    parse(p, varargin{:});

    % Retrieve the line properties
    lineColor = p.Results.lineColor;
    lineStyle = p.Results.lineStyle;
    lineWidth = p.Results.lineWidth;
    alpha = p.Results.alpha;

    % Adjust color for alpha value if necessary
    if length(lineColor) == 3  % RGB only, add alpha
        lineColor = [lineColor, alpha];
    else  % RGB and alpha
        lineColor(4) = alpha;
    end

    % Loop through each x-coordinate and draw a line
    for i = 1:length(x)
        line([x(i), x(i)], ylimit, 'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth', lineWidth);
    end
end
