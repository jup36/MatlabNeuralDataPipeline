function selectedColorC = cmapSampler(n)

% Load the autumn colormap
cmap = parula;
%cmap = flip(cmap, 1); 

% Determine the total number of colors in the colormap
numColors = size(cmap, 1);

% Exclude the first and last colors
startIdx = 2;
endIdx = numColors - 1;
effectiveColors = cmap(startIdx:endIdx, :);

% Determine the indices of four evenly spaced colors
numEffectiveColors = size(effectiveColors, 1);
indices = round(linspace(1, numEffectiveColors, n));

% Extract the four colors
selectedColors = effectiveColors(indices, :);
selectedColorC = mat2cell(selectedColors, ones(n, 1), size(selectedColors, 2)); 

% % Display the selected colors
% disp('Selected colors:');
% disp(selectedColors);
% 
% % Optionally, plot the selected colors for visualization
% figure;
% hold on;
% for i = 1:n
%     plot(1:10, i*ones(1,10), 'Color', selectedColors(i,:), 'LineWidth', 10);
% end
% hold off;
% title('Selected Colors from Autumn Colormap');
end