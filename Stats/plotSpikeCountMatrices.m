function fig = plotSpikeCountMatrices(spikeCountCell, caxislim, subplotspacing)
% caxislim = [0 5];
numCells = length(spikeCountCell); % Number of matrices/cells

% Create figure and adjust its size
fig = figure;
screenSize = get(0, 'ScreenSize'); % Get the screen size
screenHeight = screenSize(4);

figureWidth = screenHeight / 4; % Example: set figure width to a third of the screen width
figureHeight = screenHeight; % Use the full height

% Adjust figure position and size: [left, bottom, width, height]
set(fig, 'Position', [100, 100, figureWidth, figureHeight]);

% Adjust subplot spacing manually as before
% Your subplot creation and adjustment code here
ax = gobjects(numCells, 1); % Preallocate an array of graphics objects for axes
for i = 1:numCells
    ax(i) = subplot(numCells, 1, i);
    imagesc(spikeCountCell{i});
    colormap(flipud(gray));
    clim(caxislim);
    set(ax(i), 'YAxisLocation', 'left');
    if i ~= numCells
         set(ax(i), 'XTickLabel', []);
         set(gca, 'XTick', [])
    else
        set(ax(i), 'XAxisLocation', 'bottom');
        xlabel('Time Bin');
    end
    set(gca, 'TickDir', 'out')
    %ylabel(sprintf('Cell %d', i));
end

% Adjust subplot spacing
adjustSubplotSpacing(ax, numCells, subplotspacing);

%sgtitle('Spike Count Matrix Visualization');

    function adjustSubplotSpacing(ax, numCells, spacing)
        bottom = 0.1; % Adjust bottom start position
        top = 0.9; % Leave some space for title
        subplotHeight = (top - bottom) / numCells; % Calculate the height of each subplot
        %spacing = 0.002; % Adjust spacing between subplots
        for j = 1:numCells
            pos = get(ax(j), 'Position');
            pos(2) = bottom + (numCells-j)*(subplotHeight + spacing); % Adjust bottom position
            pos(4) = subplotHeight - spacing; % Adjust height to account for spacing
            set(ax(j), 'Position', pos);
        end
    end


end

