function cmap = buildPastelToVividColormap(center_color, nColor, adjustFactor)

% Define the center color in RGB
%center_color = [0.9843, 0.7059, 0.6824];
%adjustFactor = 0.75; 

% Initialize the cmap
cmap_pastel = center_color.*(1/adjustFactor); 
cmap_vivid = center_color.*adjustFactor; 

for i = 1:3
    cmap_pastel(i) = min(1, cmap_pastel(i)); 
    cmap_vivid(i) = min(1, cmap_vivid(i)); 
end

for i = 1:3 
% Calculate the step for pastel to vivid transition
    cmap(:, i) = linspace(cmap_pastel(:, i), cmap_vivid(:, i), nColor);
end

% Display the cmap
figure;
colormap(cmap);
colorbar;
