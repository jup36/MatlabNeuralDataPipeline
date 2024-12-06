function cmap = pastelColormap(nColors)
    baseColors = lines(nColors); % Start with MATLAB's 'lines' colormap
    cmap = 0.8 * baseColors + 0.2; % Blend with white to make pastel
end