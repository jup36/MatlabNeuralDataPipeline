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
