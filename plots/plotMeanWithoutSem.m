function hFig = plotMeanWithoutSem(meanData, timeX, varargin)
% varargin can provide names to be used for the figure legend

numbGroups = size(meanData, 1); 

if ~isempty(varargin)
    legendChar = sum(cell2mat(cellfun(@ischar, varargin{1}, 'UniformOutput', false))); 
    if size(meanData, 1) == legendChar
        nameLegend = false;
    else
        nameLegend = true; 
    end
else
    nameLegend = true; 
end


% Create figure and axes
hFig = figure; 
hold on;

% Get the total number of cell arrays
% (Assuming your timestamps are stored in a cell array called timeStampCells)

% Use 'cool' colormap and sample it
fullColormap = colormap('cool');
cmapIndices = round(linspace(1, size(fullColormap, 1), numbGroups));
cmap = fullColormap(cmapIndices, :);

% Define the alpha level for the shaded SEM area
alphaValue = 0.2;

% Array to store handles for the mean data plots
meanDataHandles = zeros(1, numbGroups);
legendLabel = cell(1, numbGroups); 
for i = 1:numbGroups
    if nameLegend
        legendLabel{i} = sprintf('dataSet#%d', i); 
    else
        legendLabel{i} = varargin{1}{i}; 
    end

    % Plot the mean data
    meanDataHandles(i) = plot(timeX, meanData(i, :), 'Color', cmap(i, :), 'LineWidth', 1.5);
end

xlabel('Time');
ylabel('Value');
set(gca, 'XTick', ceil(timeX(1)):1:floor(timeX(end)))
legend(meanDataHandles, legendLabel, 'Interpreter', 'none')

hold off;

end
