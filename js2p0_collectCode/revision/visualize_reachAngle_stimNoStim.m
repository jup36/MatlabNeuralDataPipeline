
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'trI')

%% prepare color
pastel2 = slanCM('Pastel2', 10); 
%plotColorListWithNumbers(pastel2); 

pastel1 = slanCM('Pastel1', 10); 
%plotColorListWithNumbers(pastel1); 

%% control trials classified as left and right trials 
controlRiI = cellfun(@(x) contains(x, 'ri'), rezCol.controlRsAlignRchAngTrId);
controlLeI = cellfun(@(x) contains(x, 'le'), rezCol.controlRsAlignRchAngTrId);

medianAngle = nanmean([rezCol.controlRsAlignRchAngRaw, rezCol.stimFullBtRchAngRaw]); % for median subtraction

rezCol.controlRsAlignRchAngRaw = rezCol.controlRsAlignRchAngRaw-medianAngle; % for median subtraction 

controlAngRi = rezCol.controlRsAlignRchAngRaw(controlRiI);
mControlAngRi = nanmean(controlAngRi); 
controlAngRi = controlAngRi(abs(controlAngRi)<70); 

controlAngLe = rezCol.controlRsAlignRchAngRaw(controlLeI); 
mControlAngLe = nanmean(controlAngLe); 
controlAngLe = controlAngLe(abs(controlAngLe)<70); 

% controlAngRi = rezCol.controlRsAlignRchAngRawMedSub(controlRiI);
% controlAngRi = controlAngRi(abs(controlAngRi)<70); 
% controlAngLe = rezCol.controlRsAlignRchAngRawMedSub(controlLeI); 
% controlAngLe = controlAngLe(abs(controlAngLe)<70); 

figRA_control = figure; % Create a new figure or get an existing figure handle
hold on; 
%plotReachAnglesOnFig(figRA_control, {controlAngRi}, {pastel2(4,:)}, {pastel2(4,:)}, 0.25);
%plotReachAnglesOnFig(figRA_control, {controlAngLe}, {pastel1(1,:)}, {pastel1(1,:)}, 0.25);
plotReachAnglesSeparateMeanOnFig(figRA_control, {controlAngRi}, {mControlAngRi}, {pastel2(4,:)}, {pastel2(4,:)}, 0.25);
plotReachAnglesSeparateMeanOnFig(figRA_control, {controlAngLe}, {mControlAngLe}, {pastel1(1,:)}, {pastel1(1,:)}, 0.25);

ylim([0 1])

print(figRA_control, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'contolLeRiRA'), '-dpdf', '-vector')

%% breakthrough trials classified as left and right trials 
stimFullBtRiI = cellfun(@(x) contains(x, 'ri'), rezCol.stimFullBtRchAngTrId);
stimFullBtLeI = cellfun(@(x) contains(x, 'le'), rezCol.stimFullBtRchAngTrId);

rezCol.stimFullBtRchAngRaw = rezCol.stimFullBtRchAngRaw-medianAngle; 

stimFullBtAngRi = rezCol.stimFullBtRchAngRaw(stimFullBtRiI);
mStimFullBtAngRi = nanmean(stimFullBtAngRi); 
stimFullBtAngRi = stimFullBtAngRi(abs(stimFullBtAngRi)<70);

stimFullBtAngLe = rezCol.stimFullBtRchAngRaw(stimFullBtLeI);
mStimFullBtAngLe = nanmean(stimFullBtAngLe); 
stimFullBtAngLe = stimFullBtAngLe(abs(stimFullBtAngLe)<70);

figRA_stim = figure; % Create a new figure or get an existing figure handle
hold on; 
%plotReachAnglesOnFig(figRA_stim, {stimFullBtAngRi}, {pastel2(4,:)}, {pastel2(4,:)}, 0.25);
%plotReachAnglesOnFig(figRA_stim, {stimFullBtAngLe}, {pastel1(1,:)}, {pastel1(1,:)}, 0.25);
plotReachAnglesSeparateMeanOnFig(figRA_stim, {stimFullBtAngRi}, {mStimFullBtAngRi}, {pastel2(4,:)}, {pastel2(4,:)}, 0.25);
plotReachAnglesSeparateMeanOnFig(figRA_stim, {stimFullBtAngLe}, {mStimFullBtAngLe}, {pastel1(1,:)}, {pastel1(1,:)}, 0.25);
ylim([0 1])

print(figRA_stim, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'fullBtLeRiRA'), '-dpdf', '-vector')

%% prepare reach angle data
% positive angle control trials
contolPosRA = rezCol.controlRsAlignRchAngRaw(rezCol.controlRsAlignRchAngRaw>0); 
contolPosRA = contolPosRA(contolPosRA<70); 
% negative angle control trials
contolNegRA = rezCol.controlRsAlignRchAngRaw(rezCol.controlRsAlignRchAngRaw<0); 
contolNegRA = contolNegRA(contolNegRA>-70); 
% absolute angle control trials
controlAbsRA = rezCol.controlRsAlignRchAng(rezCol.controlRsAlignRchAng<70); 

% positive angle stim trials
stimPosRA = rezCol.stimFullBtRchAngRaw(rezCol.stimFullBtRchAngRaw>0); 
stimPosRA = stimPosRA(stimPosRA<70); 
% negative angle stim trials
stimNegRA = rezCol.stimFullBtRchAngRaw(rezCol.stimFullBtRchAngRaw<0); 
stimNegRA = stimNegRA(stimNegRA>-70); 
% absoulte angle stim trials
stimAbsRA = rezCol.stimFullBtRchAng(rezCol.stimFullBtRchAng<70); 

%% draw reach angle vectors
figRA_control = figure; % Create a new figure or get an existing figure handle
hold on; 
plotReachAnglesSeparateMeanOnFig(figRA_control, {contolPosRA}, {pastel1(1,:)}, {pastel1(1,:)}, 0.25);
plotReachAnglesSeparateMeanOnFig(figRA_control, {contolNegRA}, {pastel2(4,:)}, {pastel2(4,:)}, 0.25);
ylim([0 1])
set(gca, 'TickDir', 'out')
%print(figRA_control, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'contolPosNegRA'), '-dpdf', '-vector')

figRA_stim = figure; 
hold on; 
plotReachAnglesSeparateMeanOnFig(figRA_stim, {stimPosRA}, {pastel1(1,:)}, {pastel1(1,:)}, 0.25);
plotReachAnglesSeparateMeanOnFig(figRA_stim, {stimNegRA}, {pastel2(4,:)}, {pastel2(4,:)}, 0.25);
ylim([0 1])
set(gca, 'TickDir', 'out')
hold off
%print(figRA_stim, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimPosNegRA'), '-dpdf', '-vector')

%% draw histogram distribution
histogramProbDense({controlAbsRA, stimAbsRA}, {pastel1(3, :), pastel2(5, :)}, 0.4)
set(gca, 'TickDir', 'out')
%set(gca, 'YScale', 'log');
xlim([0 90])

%histogramProbDensity({controlAbsRA, stimAbsRA}, {pastel1(3, :), pastel2(5, :)}, 0.4, {5, 5})
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'stimNoStimAbsRAhist'), '-dpdf', '-vector')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotReachAnglesOnFig(figObj, angleCells, colorCells, meanColorCells, lineWidth)
    % Ensure figObj is a valid figure handle
    if ~ishandle(figObj) || ~strcmp(get(figObj, 'Type'), 'figure')
        error('The first input must be a valid figure handle.');
    end
    
    % Activate the specified figure without creating a new figure or altering its properties
    set(0, 'CurrentFigure', figObj);
    
    % Validate input lengths
    if length(angleCells) ~= length(colorCells) || length(angleCells) ~= length(meanColorCells)
        error('All input cell arrays must be of the same length.');
    end
    
    % Validate if lineWidth is provided, else set to default
    if nargin < 5 || isempty(lineWidth)
        lineWidth = 1; % Default line width
    end

    % Loop through each set of angles
    for i = 1:length(angleCells)
        angles_degrees = angleCells{i}; % Get current set of angles
        color = colorCells{i}; % Get corresponding color for individual vectors
        meanColor = meanColorCells{i}; % Get corresponding color for the mean vector

        % Convert angles to radians and rotate 90 degrees counterclockwise
        angles_radians = angles_degrees * (pi / 180) + pi/2;

        % Calculate vector components for individual vectors
        x_components = cos(angles_radians);
        y_components = sin(angles_radians);

        % Plot vectors for current condition
        quiver(zeros(size(x_components)), zeros(size(y_components)), x_components, y_components, 0, 'Color', color, 'LineWidth', lineWidth);
        
        % Calculate mean angle and rotate
        mean_angle_radians = atan2(mean(sin(angles_radians)), mean(cos(angles_radians)));

        % Calculate vector component for mean vector
        mean_x_component = 1.5*cos(mean_angle_radians);
        mean_y_component = 1.5*sin(mean_angle_radians);
        
        % Plot mean vector with twice the line width and without arrowhead using 'plot'
        plot([0, mean_x_component], [0, mean_y_component], 'Color', meanColor, 'LineWidth', lineWidth * 20);
    end
end
