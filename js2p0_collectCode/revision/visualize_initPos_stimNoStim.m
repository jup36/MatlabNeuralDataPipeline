
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI')

%% prepare color
pastel2 = slanCM('Pastel2', 10); 
%plotColorListWithNumbers(pastel2); 

pastel1 = slanCM('Pastel1', 10); 
%plotColorListWithNumbers(pastel1); 

color_L_seed = pastel1(1,:); 
color_R_seed = pastel2(4,:); 

% cmap = colormap('cool'); 
% cyan = cmap(1, :); 
% purple = cmap(end, :); 

%% generate kernel estimate distribution plots
h_noStimL = kde_2d_distribution_customColor(cell2mat(rezCol.trjNoStimRsAlignInitXY_LTrial), [1 1 1], color_L_seed, 50); 
set(gca, 'TickDir', 'out')
print(h_noStimL, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_noStimRsAlign_LTrial'), '-dpdf', '-vector')

h_noStimR = kde_2d_distribution_customColor(cell2mat(rezCol.trjNoStimRsAlignInitXY_RTrial), [1 1 1], color_R_seed, 50); 
set(gca, 'TickDir', 'out')
print(h_noStimR, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_noStimRsAlign_RTrial'), '-dpdf', '-vector')

h_stimL = kde_2d_distribution_customColor(cell2mat(rezCol.trjStimFullBtInitXY_LTrial), [1 1 1], color_L_seed, 50); 
set(gca, 'TickDir', 'out')
print(h_stimL, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_stimFullBt_LTrial'), '-dpdf', '-vector')

h_stimR = kde_2d_distribution_customColor(cell2mat(rezCol.trjStimFullBtInitXY_RTrial), [1 1 1], color_R_seed, 50); 
set(gca, 'TickDir', 'out')
print(h_stimR, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/InitPosXY_stimFullBt_RTrial'), '-dpdf', '-vector')

%% Stat
trjNoStimRsAlignInitXY_LTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_LTrial); 
trjNoStimRsAlignInitXY_RTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_RTrial); 
[~, stat.trjNoStimRsAlignInitX.p, ~, stat.trjNoStimRsAlignInitX.stats] = ttest2(trjNoStimRsAlignInitXY_LTrial(1, :), trjNoStimRsAlignInitXY_RTrial(1, :)); 
[~, stat.trjNoStimRsAlignInitY.p, ~, stat.trjNoStimRsAlignInitY.stats] = ttest2(trjNoStimRsAlignInitXY_LTrial(2, :), trjNoStimRsAlignInitXY_RTrial(2, :)); 

trjStimFullBtInitXY_LTrial = cell2mat(rezCol.trjStimFullBtInitXY_LTrial); 
trjStimFullBtInitXY_RTrial = cell2mat(rezCol.trjStimFullBtInitXY_RTrial); 
[~, stat.trjStimFullBtInitXYInitX.p, ~, stat.trjStimFullBtInitX.stats] = ttest2(trjStimFullBtInitXY_LTrial(1, :), trjStimFullBtInitXY_RTrial(1, :)); 
[~, stat.trjStimFullBtInitXYInitY.p, ~, stat.trjStimFullBtInitY.stats] = ttest2(trjStimFullBtInitXY_LTrial(2, :), trjStimFullBtInitXY_RTrial(2, :)); 

% compare breakthrough and control trials
% lateral position le and ri trials combined
[~, stat.controlVsFullBtInitX.p, ~, stat.controlVsFullBtInitX.stats] = ttest2([trjNoStimRsAlignInitXY_LTrial(1, :), trjNoStimRsAlignInitXY_RTrial(1, :)], ...
    [trjStimFullBtInitXY_LTrial(1, :), trjStimFullBtInitXY_RTrial(1, :)]); 
% lateral position le only 
[~, stat.controlVsFullBtInitX_le.p, ~, stat.controlVsFullBtInitX_le.stats] = ttest2(trjNoStimRsAlignInitXY_LTrial(1, :), ...
    trjStimFullBtInitXY_LTrial(1, :)); 
% lateral position ri only 
[~, stat.controlVsFullBtInitX_ri.p, ~, stat.controlVsFullBtInitX_ri.stats] = ttest2(trjNoStimRsAlignInitXY_RTrial(1, :), ...
    trjStimFullBtInitXY_RTrial(1, :)); 

% outward position le and ri trials combined
[~, stat.controlVsFullBtInitY.p, ~, stat.controlVsFullBtInitY.stats] = ttest2([trjNoStimRsAlignInitXY_LTrial(2, :), trjNoStimRsAlignInitXY_RTrial(2, :)], ...
    [trjStimFullBtInitXY_LTrial(2, :), trjStimFullBtInitXY_RTrial(2, :)]); 
% lateral position le only 
[~, stat.controlVsFullBtInitY_le.p, ~, stat.controlVsFullBtInitY_le.stats] = ttest2(trjNoStimRsAlignInitXY_LTrial(2, :), ...
    trjStimFullBtInitXY_LTrial(2, :)); 
% lateral position ri only 
[~, stat.controlVsFullBtInitY_ri.p, ~, stat.controlVsFullBtInitY_ri.stats] = ttest2(trjNoStimRsAlignInitXY_RTrial(2, :), ...
    trjStimFullBtInitXY_RTrial(2, :)); 



%% Calculate symmetric initial position coverage
trjNoStimRsAlignInitXY_LTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_LTrial); 
trjNoStimRsAlignInitXY_RTrial = cell2mat(rezCol.trjNoStimRsAlignInitXY_RTrial); 

trjNoStimRsAlignInitX_LTrial = trjNoStimRsAlignInitXY_LTrial(1, :); 
trjNoStimRsAlignInitX_RTrial = trjNoStimRsAlignInitXY_RTrial(1, :);

rangesX = calculateSymmetricDataCoverage([trjNoStimRsAlignInitX_LTrial, trjNoStimRsAlignInitX_RTrial], 0.1); 

trjNoStimRsAlignInitY_LTrial = trjNoStimRsAlignInitXY_LTrial(2, :); 
trjNoStimRsAlignInitY_RTrial = trjNoStimRsAlignInitXY_RTrial(2, :);

rangesY = calculateSymmetricDataCoverage([trjNoStimRsAlignInitY_LTrial, trjNoStimRsAlignInitY_RTrial], 0.1); 

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

function ranges = calculatePercentileRanges(data, step)
    % Calculate percentile ranges for given step sizes in a dataset.
    %
    % Parameters:
    %   data - An array containing the data points.
    %   step - The step size for the percentiles (e.g., 10 for 10%, 20% ...).
    %
    % Returns:
    %   ranges - A matrix where each row represents the [lower_bound, upper_bound] of the percentile range.

    if nargin < 2
        step = 10; % Default step size is 10 if not specified
    end
    
    percentiles = 0:step:100; % Create a vector of percentiles from 0 to 100
    values = prctile(data, percentiles); % Calculate the percentiles
    
    % Initialize the ranges matrix
    numRanges = length(values) - 1;
    ranges = zeros(numRanges, 2);
    
    for i = 1:numRanges
        ranges(i, 1) = values(i);
        ranges(i, 2) = values(i + 1);
    end
end
