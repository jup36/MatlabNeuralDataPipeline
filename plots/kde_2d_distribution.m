function h = kde_2d_distribution(data)

% Assuming data is your 2-by-N matrix where the first row is x and the second row is y
x = data(1,:);
y = data(2,:);

% Define grid for evaluation
[X, Y] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));

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


end