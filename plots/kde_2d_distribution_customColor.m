function h = kde_2d_distribution_customColor(data, seedColor, targetColor, N)
% Visualizes the 2D kernel density estimate of the provided data with a custom colormap.
%
% Parameters:
% data - A 2-by-N matrix of data points (first row for x values, second row for y values).
% seedColor - The RGB color to start the colormap with.
% targetColor - The RGB color to end the colormap with.
% N - The number of steps in the colormap.

% Validate input colors
if size(seedColor, 2) ~= 3 || size(targetColor, 2) ~= 3
    error('Seed and target colors must be 1x3 RGB vectors.');
end

% Assuming data is your 2-by-N matrix where the first row is x and the second row is y
x = data(1,:);
y = data(2,:);

% Define grid for evaluation
%[X, Y] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
[X, Y] = meshgrid(linspace(-10, 10, 100), linspace(-10, 10, 100));

% Flatten X and Y for ksdensity input
gridX = X(:);
gridY = Y(:);

% Perform kernel density estimation
[f, ~] = ksdensity([x', y'], [gridX gridY]);

% Reshape the output to fit the meshgrid dimensions
F = reshape(f, size(X));

% Plotting the density
h = figure;
contourf(X, Y, F, 50, 'LineColor', 'none'); % Adjust number of contour levels as needed
colorbar; % Shows the color scale
title('Density Distribution of 2D Points');
xlabel('X');
ylabel('Y');
xlim([-5 5])
ylim([-5 5])
clim([0 0.05])

% Generate and apply the custom colormap
cmap = customColormap(seedColor, targetColor, N);
colormap(cmap);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = customColormap(seed, target, N)
    % customColormap Generates a custom colormap.
    %
    %   cmap = customColormap(seed, target, N) creates a colormap with N colors,
    %   starting from 'seed' color and transitioning to 'target' color.
    %   'seed' and 'target' are RGB colors specified as 1x3 vectors with
    %   values in the range [0, 1]. 'N' is the total number of colors in the
    %   colormap.

    % Initialize the colormap matrix
    cmap = zeros(N, 3);
    
    % Linearly interpolate between seed and target colors
    for i = 1:3 % For R, G, and B components
        cmap(:, i) = linspace(seed(i), target(i), N);
    end
end


end
