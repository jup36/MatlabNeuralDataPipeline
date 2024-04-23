
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI')

%% prepare color
pastel2 = slanCM('Pastel2', 10); 
plotColorListWithNumbers(pastel2); 

pastel1 = slanCM('Pastel1', 10); 
plotColorListWithNumbers(pastel1); 

color_L_seed = pastel1(1,:); 
color_R_seed = pastel2(4,:); 

% cmap = colormap('cool'); 
% cyan = cmap(1, :); 
% purple = cmap(end, :); 

%% generate kernel estimate distribution plots
h_noStimL = kde_2d_distribution_customColor(cell2mat(rezCol.trjNoStimRsAlignInitXY_LTrial), [1 1 1], color_L, 50); 
set(gca, 'TickDir', 'out')
print(h_noStimL, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_noStimRsAlign_LTrial'), '-dpdf', '-vector')

h_noStimR = kde_2d_distribution_customColor(cell2mat(rezCol.trjNoStimRsAlignInitXY_RTrial), [1 1 1], color_R, 50); 
set(gca, 'TickDir', 'out')
print(h_noStimR, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_noStimRsAlign_RTrial'), '-dpdf', '-vector')

h_stimL = kde_2d_distribution_customColor(cell2mat(rezCol.trjStimFullBtInitXY_LTrial), [1 1 1], color_L, 50); 
set(gca, 'TickDir', 'out')
print(h_stimL, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_stimFullBt_LTrial'), '-dpdf', '-vector')

h_stimR = kde_2d_distribution_customColor(cell2mat(rezCol.trjStimFullBtInitXY_RTrial), [1 1 1], color_R, 50); 
set(gca, 'TickDir', 'out')
print(h_stimR, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_stimFullBt_RTrial'), '-dpdf', '-vector')

%% Stat
trjNoStimRsAlignInitXY_LTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_LTrial); 
trjNoStimRsAlignInitXY_RTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_RTrial); 
[~, stat.trjNoStimRsAlignInitX.p, ~, stat.trjNoStimRsAlignInitX.stats] = ttest2(trjNoStimRsAlignInitXY_LTrial(1, :), trjNoStimRsAlignInitXY_RTrial(1, :)); 

trjStimFullBtInitXY_LTrial = cell2mat(rezCol.trjStimFullBtInitXY_LTrial); 
trjStimFullBtInitXY_RTrial = cell2mat(rezCol.trjStimFullBtInitXY_RTrial); 
[~, stat.trjStimFullBtInitXYInitX.p, ~, stat.trjStimFullBtInitX.stats] = ttest2(trjStimFullBtInitXY_LTrial(1, :), trjStimFullBtInitXY_RTrial(1, :)); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = kde_2d_distribution(data, customColorMap)

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