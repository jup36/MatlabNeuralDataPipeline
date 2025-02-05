
filePaths = {'/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_121324/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_121624/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_121724/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_121824/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_121924/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_122024/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_122324/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_122424/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_122624/task', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda/m1049_122724/task'
             };

figSaveDir = fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda', 'collectFigure'); 
if exist(figSaveDir, "dir")~=7
    mkdir(figSaveDir)
end

mID = regexp(filePaths{1}, 'm\d{1,4}', 'match', 'once');

% prepare colormaps
blueShades = generateColormap([176 226 255]./255, [0 0 128]./255, 200); % for blocks
cool = colormap('cool'); % for Go/No-Go
pastels = slanCM('Pastel1', 7); % for sigD data

dPrmC = cell(length(filePaths), 2); 
dPrmTot = zeros(length(filePaths), 1); 

%% Main Loop
for f = 1:length(filePaths) 
    % load rez of LickAnalysis
    %[~, header] = fileparts(filePaths{f});
    header = regexp(filePaths{f}, 'm\d{1,4}_\d{6}', 'match', 'once'); 
    load(fullfile(filePaths{f}, 'Matfiles', [header, '_LickAnalysis']), 'rez', 'var');
    
    dPrmTot(f, 1) = rez.sigD.dprime; 

    % collect dPrime data 
    dPrmC{f, 1} = header; 
    dPrmC{f, 2} = cell2mat(cellfun(@(a) a.dprime, rez.sigDBlocks, 'UniformOutput', false)); 

    % collect latency data 
    latGoC{f, 1} = header; 
    latGoC{f, 2} = rez.lat.rwdFstLatBlockMean; 

    latNogoC{f, 1} = header;
    latNogoC{f, 2} = rez.lat.pnsFstLatBlockMean; 
    
    fprintf('Completed file #%d\n', f); 
end

%% plot
% plot dPrime
hfig = sigDrezDprmPlotAcrossSession('dPrime Across Session', dPrmC(:, 2), blueShades); 
print(hfig, fullfile(figSaveDir, [mID, '_LickAnalysis_acrossSession_dPrime']), '-dpdf', '-vector');

% plot first lick latency
hLat = firstLatencyPlotGngAcrossSession('first lick latency across session', latGoC(:, 2), latNogoC(:, 2), cool); 
print(hLat, fullfile(figSaveDir, [mID, '_LickAnalysis_acrossSession_firstLickLatencyGoNogo']), '-dpdf', '-vector');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% z-score normalization
function [zrez, zvar] = binnedLickZWrapper(var, lickTimeC, varargin)
%This function takes lick timestamps across trials in a cell array. The timestamps are
% binned, counted, and z-score normalized. Trial-group variables can be
% provided as varargin as a cell array containing trial numbers or
% logicals.

% lickTimes with trial grouping
if ~isempty(varargin) % use grouping
    zvar.trGroups = varargin{1};
    lickTimes = cellfun(@(a) lickTimeC(a), zvar.trGroups, 'UniformOutput', false);
else % no grouping
    zvar.trGroups = {true(length(lickTimeC), 1)}; % there's only one group
    lickTimes = lickTimeC;
end

% observed lick rasters
%rasterPlotCell(lickTimes, 0.5);

%% normalize across all trials
lickTimesConcat = [lickTimes{:}];
lickTimesSort = cell(1, length(lickTimesConcat));
for g = 1:length(zvar.trGroups)
    if islogical(zvar.trGroups{g})
        trGroupsI = find(zvar.trGroups{g});
    else
        trGroupsI = zvar.trGroups{g};
    end
    lickTimesSort(trGroupsI) = deal(lickTimes{g});
end

[zrez.lickZ, zvar.binEdgesLickZ] = binTimestampsZscore(lickTimesSort, var);
zrez.smLickZ = smooth2a(zrez.lickZ, 0, 3); % smooth across time

zvar.timeX_binnedLickZ = zvar.binEdgesLickZ(1:end-1)+ var.binSize/2;

%% normalize within each group
zrez.lickGroupZ = binTimestampsGroupZscore(lickTimes, var);
zrez.smLickGroupZ = cellfun(@(a) smooth2a(a, 0, 3), zrez.lickGroupZ, 'UniformOutput', false);

%% mean, sem of lickZ, lickGroupZ
zrez.mLickZ = zeros(length(zvar.trGroups), size(zrez.lickZ, 2));
zrez.sLickZ = zeros(length(zvar.trGroups), size(zrez.lickZ, 2));

zrez.mLickGroupZ = zeros(length(zvar.trGroups), size(zrez.lickZ, 2));
zrez.sLickGroupZ = zeros(length(zvar.trGroups), size(zrez.lickZ, 2));

for group = 1:length(zvar.trGroups)
    [zrez.mLickZ(group, :), ~, zrez.sLickZ(group, :)] = meanstdsem(zrez.smLickZ(zvar.trGroups{group}(:), :));
    [zrez.mLickGroupZ(group, :), ~, zrez.sLickGroupZ(group, :)] = meanstdsem(zrez.smLickGroupZ{group});
end


% plotMeanSem(zrez.mLickZ, zrez.sLickZ, zvar.timeX_binnedLickZ) % plot the across-trial z-score licks
% plotMeanSem(zrez.mLickGroupZ, zrez.sLickGroupZ, zvar.timeX_binnedLickZ) % plot the across-trial z-score licks


    function zscoreC = binTimestampsGroupZscore(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numGroups = length(timeStampC);

        for gr = 1:numGroups
            zscoreC{gr} = binTimestampsZscore(timeStampC{gr}, var);
        end

        %zscoreMat = cell2mat([zscoreC(:)]);

    end


    function [z_score_binCounts_trials, binEdges] = binTimestampsZscore(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numTrials = length(timeStampC);

        % Define bin edges
        binEdges = var.minMaxTpre(1):var.binSize:var.minMaxTpost(2); % e.g. from -3s to 8s in 100ms (0.1s) increments

        baselineLogic = binEdges < 0;

        % Initialize storage for binned and z-score normalized lick counts
        binCountsTrials = zeros(length(timeStampC), length(binEdges) - 1);

        for i = 1:numTrials
            % Bin the lick timestamps for the trial
            if ~isempty(timeStampC{i})
                binCounts = histc(timeStampC{i}, binEdges);

                % Exclude the last bin
                binCounts = binCounts(1:end-1);
                binCountsTrials(i, :) = binCounts;
            end
        end

        % Calculate mean and standard deviation from the baseline period across all trials (first 50 bins)
        baseline_mean = mean(mean(binCountsTrials(:, baselineLogic)));
        baseline_std = std(binCountsTrials(:, baselineLogic), [], 'all'); % the 'all' option computes the std considering all elements

        % Z-score normalization across all trials
        z_score_binCounts_trials = (binCountsTrials - baseline_mean) / baseline_std;

    end
end


%% Signal Detection Theoretic analysis of licking
function drez = signalDetectionLickAnalysis(tbytDat, trRwdI, trPnsI)
%trRwdI = var.stimopts.rewarded_stim; %
%trPnsI = var.stimopts.punished_stim;

% classify trials
rewarded = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.water}, 'un', 0));
punished = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.airpuff}, 'un', 0));

% compute overall hit, miss, fa, cr rates
rawHitRate = sum(rewarded(:) & trRwdI(:))/sum(trRwdI); % raw hit rate
drez.hitRate = adjustRate(rawHitRate, sum(trRwdI));

drez.missRate = 1 - drez.hitRate; % miss rate

if sum(trPnsI) > 5 % If there were No-Go trials
    rawFaRate = sum(punished(:) & trPnsI(:))/sum(trPnsI); % raw false alarm rate
    drez.FaRate = adjustRate(rawFaRate, sum(trPnsI));

    drez.CrRate = 1 - drez.FaRate; % correct rejection rate

    % d' (discriminability): The metric d' is a measure of how well an
    % observer can distinguish between two different stimuli or
    % conditions (signal vs. noise).
    drez.dprime = norminv(drez.hitRate) - norminv(drez.FaRate);
end

    function adjustedRate = adjustRate(rawRate, N)
        % Adjusts hit or false alarm rate if it's 1 or 0
        if rawRate == 1
            adjustedRate = (N-0.5) / N;
        elseif rawRate == 0
            adjustedRate = 0.5 / N;
        else
            adjustedRate = rawRate;
        end
    end

end



function [lat] = latencyOfFirstLicks(tbytDat, trRwdI, trPnsI, cueWin)

% First lick latency across blocks
blocks = divideTrials(length(tbytDat), 3, 10, 10);

lickTimeC = cellfun(@(a, b) a-b, {tbytDat.Lick}, {tbytDat.evtOn}, 'un', 0);

lickTimeCueC = cellfun(@(a) a(a>=cueWin(1) & a<=cueWin(2)), lickTimeC, 'UniformOutput', false);

lickCounts = cell2mat(cellfun(@length, lickTimeCueC, 'UniformOutput', false));
valTrI = lickCounts > 0;
valRwdTrI = valTrI(:) & trRwdI(:);

% First lick latency
lat.allFstLat = cell(1, length(tbytDat)); lat.allFstLat(:) = {NaN};
lat.rwdFstLat = cell(1, length(tbytDat)); lat.rwdFstLat(:) = {NaN};
lat.allFstLat(valTrI) = cellfun(@(a) a(1), lickTimeCueC(valTrI), 'UniformOutput', false);
lat.rwdFstLat(valRwdTrI) = cellfun(@(a) a(1), lickTimeCueC(valRwdTrI), 'UniformOutput', false);

% across blocks
allFstLatBlocks = cellfun(@(a) cell2mat(lat.allFstLat(a)), blocks, 'UniformOutput', false);
lat.allFstLatBlockMean = cell2mat(cellfun(@nanmean, allFstLatBlocks, 'UniformOutput', false));

rwdFstLatBlocks = cellfun(@(a) cell2mat(lat.rwdFstLat(a)), blocks, 'UniformOutput', false);
lat.rwdFstLatBlockMean = cell2mat(cellfun(@nanmean, rwdFstLatBlocks, 'UniformOutput', false));

% for No-Go trials
if sum(trPnsI)>5 % if there are No-Go trials
    valPnsTrI = valTrI(:) & trPnsI(:);
    lat.pnsFstLat = cell(1, length(tbytDat)); lat.pnsFstLat(:) = {NaN};
    lat.pnsFstLat(valPnsTrI) = cellfun(@(a) a(1), lickTimeCueC(valPnsTrI), 'UniformOutput', false);
    pnsFstLatBlocks = cellfun(@(a) cell2mat(lat.pnsFstLat(a)), blocks, 'UniformOutput', false);
    lat.pnsFstLatBlockMean = cell2mat(cellfun(@nanmean, pnsFstLatBlocks, 'UniformOutput', false));
end

end


%% examine and compare the simulated versus observed lick data
function rez = examineLickCWrapper(lickSimInC, lickTimeSimC, minMaxTpre, minMaxTcue, minMaxTpost)
for j = 1:size(lickSimInC, 1) % increment blocks
    % pre-cue
    [rez.lickCntAct(j, 1), rez.lickIntAct(j, 1), rez.lickLambdaAct(j, 1), rez.lickIntActWhole{j, 1}] = examineLickC(lickSimInC{j, 1}, minMaxTpre);
    [rez.lickCntSim(j, 1), rez.lickIntSim(j, 1), rez.lickLambdaSim(j, 1), rez.lickIntSimWhole{j, 1}] = examineLickC(lickTimeSimC{j, 1}, minMaxTpre);

    % cue
    [rez.lickCntAct(j, 2), rez.lickIntAct(j, 2), rez.lickLambdaAct(j, 2), rez.lickIntActWhole{j, 2}] = examineLickC(lickSimInC{j, 2}, minMaxTcue);
    [rez.lickCntSim(j, 2), rez.lickIntSim(j, 2), rez.lickLambdaSim(j, 2), rez.lickIntSimWhole{j, 2}] = examineLickC(lickTimeSimC{j, 2}, minMaxTcue);

    % post-cue
    [rez.lickCntAct(j, 3), rez.lickIntAct(j, 3), rez.lickLambdaAct(j, 3), rez.lickIntActWhole{j, 3}] = examineLickC(lickSimInC{j, 3}, minMaxTpost);
    [rez.lickCntSim(j, 3), rez.lickIntSim(j, 3), rez.lickLambdaSim(j, 3), rez.lickIntSimWhole{j, 3}] = examineLickC(lickTimeSimC{j, 3}, minMaxTpost);
end

    function [lickCnt, lickInt, lickLambda] = examineLickC(lickC, minMaxTime)

        lickCnt = cellfun(@(a) length(a)/abs(diff(minMaxTime)), lickC, 'UniformOutput', false);
        lickInt = nanmean(cell2mat(cellfun(@(a) diff([minMaxTime(1); a]), lickC, 'UniformOutput', false)'));
        lickLambda = 1/lickInt;

    end
end


function hfig = rasterPlotLickSim(lickSimCell, numTrs)
%lickSimCell contains lick timestamps sampled from the exponential
% distribution fitted to the data

numSamples = unique(cell2mat(cellfun(@length, lickSimCell, 'UniformOutput', false)));
trs = randperm(numSamples, numTrs);

lickC = cell(size(lickSimCell, 1), 1); % blocks
for ts = 1:size(lickSimCell, 1)
    % choose block (early, middle, late)
    tsLicksC = cellfun(@(a) a(trs), lickSimCell(ts, :), 'UniformOutput', false);

    % concatenate epochs (pre-cue, cue, post-cue)
    tsLickTsC = cellfun(@(a, b, c) [a; b; c], tsLicksC{1}, tsLicksC{2}, tsLicksC{3}, 'UniformOutput', false); % concatenate epochs
    lickC{ts} = tsLickTsC(:);
end

hfig = rasterPlotCell(lickC);


    function hfig = rasterPlotCell(cellArrays)
        % RASTERPLOTMULTI Generates a raster plot for multiple cell arrays
        % Each entry in the cellArrays is a cell array containing cells representing trials with timestamps.

        % Get the total number of cell arrays
        numCellArrays = length(cellArrays);

        % Use 'cool' colormap and sample it
        fullColormap = colormap('cool');
        cmapIndices = round(linspace(1, size(fullColormap, 1), numCellArrays));
        cmap = fullColormap(cmapIndices, :);

        % Prepare a figure for the raster plot
        hfig = figure;
        hold on;

        currentTrialOffset = 0;

        % Loop through each cell array
        for c = 1:numCellArrays
            cellArray = cellArrays{c};
            numTrials = length(cellArray);
            currentColor = cmap(c, :); % Get the color for the current cell array

            % Loop through each trial and plot the rasters
            for trial = 1:numTrials
                % Retrieve the time points for the current trial
                timePoints = cellArray{trial};

                % Plot a line for each time point in the current trial using the line function
                for tp = 1:length(timePoints)
                    line([timePoints(tp), timePoints(tp)], [trial + currentTrialOffset - 0.7, trial + currentTrialOffset + 0.7], 'Color', currentColor, 'LineWidth', 2);
                end
            end

            % Update the current trial offset
            currentTrialOffset = currentTrialOffset + numTrials;
        end

        % Set the y-axis to display trials from top to bottom
        set(gca, 'YDir', 'reverse');
        xlim([floor(min(cell2mat(cellfun(@(x) min(cell2mat(x)), cellArrays, 'UniformOutput', false)))), ceil(max(cell2mat(cellfun(@(x) max(cell2mat(x)), cellArrays, 'UniformOutput', false))))]);
        ylim([0, currentTrialOffset + 1]);
        xlabel('Time');
        ylabel('Trial #');
        title('Raster Plot');
        hold off;
    end
end


function timeStampOut = lickSim(timeStampC, minMaxTime, numReps)

% timeStampC: a cell array that contains timestamps of which interval
% times to be modeled.
% timeStampC = lickTimeEarlyCue;
% minMaxTime = [0 5];
% minTime = 0; % sec
% numReps = 1000;

% sanity check: ensure that all timestamps are within the min and max range
timeStampC = cellfun(@(a) a(a >= minMaxTime(1) & a <= minMaxTime(2)), timeStampC, 'UniformOutput', false);

intTimeC = cellfun(@(a) diff([minTime; a]), timeStampC, 'UniformOutput', false); % transform to interval time
intTime = cell2mat(intTimeC'); % take all interval times

% fit exponential distribution
expd = fitdist(intTime, 'Exponential');

lambda = 1/mean(intTime); % lambda is the rate parameter of the exponential distribution

timeStampOut = cell(1, numReps);
% generate timestamps
for jj = 1:numReps
    timeStampOut{jj} = expTimestampsTrial(expd, minMaxTime);
end



    function timeStamps = expTimestampsTrial(expd, minMaxTime)
        minTime = min(minMaxTime);
        maxTime = max(minMaxTime);

        curTime = minTime; % starting point
        intvs = []; % generated intervals
        while curTime <= maxTime
            intv = random(expd, 1, 1);
            curTime = curTime + intv;
            if curTime <= maxTime
                intvs = [intvs; intv];
            end
        end
        % generate cumulative timestamps
        timeStamps = cumsum(intvs) + minTime;
    end

end


function rasterPlot(cellArray, colorArray, alphaValue)
% RASTERPLOT Generates a raster plot for the given cell array
% Each cell in the cellArray represents a trial and contains time points
% relative to the trial onset.
% colorArray is an RGB array for setting the color of the rasters.

% Determine the number of trials
numTrials = length(cellArray);

% Check if colorArray length matches cellArray length
if size(colorArray, 1) == 1
    colorArray = repmat(colorArray, [numTrials, 1]);
elseif size(colorArray, 1) ~= numTrials
    error('The length of the color array should match the number of trials.');
end

% Prepare a figure for the raster plot
figure;
hold on;

% Loop through each trial and plot the rasters
for trial = 1:numTrials
    % Retrieve the time points for the current trial
    timePoints = cellArray{trial};
    currentColor = [colorArray(trial, :) alphaValue]; % Get the color for the current trial and append alpha value

    % Plot a line for each time point in the current trial using the line function
    for tp = 1:length(timePoints)
        line([timePoints(tp), timePoints(tp)], [trial-0.7, trial+0.7], 'Color', currentColor, 'LineWidth', 2);
    end
end

% Set the y-axis to display trials from top to bottom
set(gca, 'YDir', 'reverse');
xlim([floor(min(cell2mat(cellfun(@min, cellArray, 'UniformOutput', false)'))), ceil(max(cell2mat(cellfun(@max, cellArray, 'UniformOutput', false)')))]);
ylim([0, numTrials + 1]);
xlabel('Time');
ylabel('Trial #');
title('Raster Plot');
hold off;
end


function rasterPlotCell(cellArrays)
% RASTERPLOTMULTI Generates a raster plot for multiple cell arrays
% Each entry in the cellArrays is a cell array containing cells representing trials with timestamps.

% Get the total number of cell arrays
numCellArrays = length(cellArrays);

% Use 'cool' colormap and sample it
fullColormap = colormap('cool');
cmapIndices = round(linspace(1, size(fullColormap, 1), numCellArrays));
cmap = fullColormap(cmapIndices, :);

% Prepare a figure for the raster plot
figure;
hold on;

currentTrialOffset = 0;

% Loop through each cell array
for c = 1:numCellArrays
    cellArray = cellArrays{c};
    numTrials = length(cellArray);
    currentColor = cmap(c, :); % Get the color for the current cell array

    % Loop through each trial and plot the rasters
    for trial = 1:numTrials
        % Retrieve the time points for the current trial
        timePoints = cellArray{trial};

        % Plot a line for each time point in the current trial using the line function
        for tp = 1:length(timePoints)
            line([timePoints(tp), timePoints(tp)], [trial + currentTrialOffset - 0.7, trial + currentTrialOffset + 0.7], 'Color', currentColor, 'LineWidth', 2);
        end
    end

    % Update the current trial offset
    currentTrialOffset = currentTrialOffset + numTrials;
end

% Set the y-axis to display trials from top to bottom
set(gca, 'YDir', 'reverse');
xlim([floor(min(cell2mat(cellfun(@(x) min(cell2mat(x)), cellArrays, 'UniformOutput', false)))), ceil(max(cell2mat(cellfun(@(x) max(cell2mat(x)), cellArrays, 'UniformOutput', false))))]);
ylim([0, currentTrialOffset + 1]);
xlabel('Time');
ylabel('Trial #');
title('Raster Plot');
hold off;
end


function plotIntervalTimes(intTime, expd)
% intTime: An array comprising interval times
minInt = floor(min(intTime));
maxInt = ceil(max(intTime));

figure;
histogram(intTime, 'Normalization', 'pdf', 'DisplayStyle', 'stairs'); % original data
hold on;
x = linspace(minInt, maxInt, 100);
y = pdf(expd, x);
plot(x, y, 'r', 'LineWidth', 2); % fitted distribution
xlabel('Interval Times');
ylabel('Probability Density');
legend('Original Data', 'Fitted Distribution');
hold off;
end


function h = blockMeanPlot(blockMeanC)
% blockMeanC = {lat.rwdFstLatBlockMean; lat.pnsFstLatBlockMean};
% Use 'cool' colormap and sample it
fullColormap = colormap('cool');
cmapIndices = round(linspace(1, size(fullColormap, 1), length(blockMeanC)));
cmap = fullColormap(cmapIndices, :);

h = figure; hold on;
for i = 1:length(blockMeanC)
    scatter(1:length(blockMeanC{i}), blockMeanC{i}, 100, cmap(i, :), 'filled')
    plot(1:length(blockMeanC{i}), blockMeanC{i}, 'Color', cmap(i, :), 'LineStyle', ':')
end
% Determine the full range of x-values and y-values
allLengths = cellfun(@length, blockMeanC);
allValues = cell2mat(blockMeanC);

minX = 1;
maxX = max(allLengths);
minY = min(allValues(:));
maxY = max(allValues(:));

% Compute 10% margins for x and y axes
xMargin = 0.1 * (maxX - minX);
yMargin = 0.1 * (maxY - minY);

% Update xlim and ylim with the new ranges
xlim([minX - xMargin, maxX + xMargin]);
ylim([minY - yMargin, maxY + yMargin]);

% Set TickDir to 'out', YTick to have six ticks, and YTickLabels with one decimal point
xTicks = 1:length(blockMeanC{i});
yTicks = linspace(minY - yMargin, maxY + yMargin, 6);
yTickLabels = arrayfun(@(x) sprintf('%.1f', x), yTicks, 'UniformOutput', false);

set(gca, 'TickDir', 'out', 'XTick', xTicks, 'YTick', yTicks, 'YTickLabel', yTickLabels, 'XLabel', 'Blocks', 'YLabel', 'Latency for first lick (s)');
hold off;
%print(h, fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101723', 'Figure', 'firstLickLatencyBlocks'), '-dpdf', '-vector');  % '-painters' ensures the output is vector graphics
end


function [lickSimInC, lickSimOutC, lickCntC] = lickSimWrapper(var, rez)
% allocate
lickSimInC = cell(length(var.trSets), 3); % #trSets-by-#epochs, each cell to contain timestamps of # reps
lickSimOutC = cell(length(var.trSets), 3, 2); % #trSets-by-#epochs-by-2, each cell to contain timestamps of # reps
lickCntC = cell(length(var.trSets), 3); % #trSets-by-#epochs, each cell to contain

% divide epochs (pre-, post-, and cue epochs)
for ts = 1:size(lickSimOutC, 1) % increment trial sets (blocks)
    % pre-cue epoch
    lickSimInC{ts, 1} = cellfun(@(a) a(a < 0), rez.lickTimeC(var.trSets{ts}), 'UniformOutput', false); % pre-cue epoch of the current trial set
    lickCntC{ts, 1} = cellfun(@length, lickSimInC{ts, 1}, 'UniformOutput', false);
    [lickSimOutC{ts, 1, 1}, lickSimOutC{ts, 1, 2}] = lickSimPoissonExpBin(lickSimInC{ts, 1}, var.minMaxTpre, var.simBinSize, var.simRep);

    % cue epoch
    lickSimInC{ts, 2} = cellfun(@(a) a(a >= 0 & a <= var.trDur), rez.lickTimeC(var.trSets{ts}), 'UniformOutput', false); % cue epoch of the current trial set
    lickCntC{ts, 2} = cellfun(@length, lickSimInC{ts, 2}, 'UniformOutput', false);
    [lickSimOutC{ts, 2, 1}, lickSimOutC{ts, 2, 2}] = lickSimPoissonExpBin(lickSimInC{ts, 2}, var.minMaxTcue, var.simBinSize, var.simRep);

    % post-cue epoch
    lickSimInC{ts, 3} = cellfun(@(a) a(a > var.trDur), rez.lickTimeC(var.trSets{ts}), 'UniformOutput', false); % cue epoch of the current trial set
    lickCntC{ts, 3} = cellfun(@length, lickSimInC{ts, 3}, 'UniformOutput', false);
    [lickSimOutC{ts, 3, 1}, lickSimOutC{ts, 3, 2}] = lickSimPoissonExpBin(lickSimInC{ts, 3}, var.minMaxTpost, var.simBinSize, var.simRep);
end


    function [timeStampOut, poisd, expd] = lickSimPoissonExpBin(timeStampC, minMaxTime, binSize, numReps)
        % Ensure timestamps are within the specified range
        timeStampC = cellfun(@(a) a(a >= minMaxTime(1) & a <= minMaxTime(2)), timeStampC, 'UniformOutput', false);

        % Compute intervals and fit exponential pdf
        intTimeC = cellfun(@(a) diff([minMaxTime(1); a]), timeStampC, 'UniformOutput', false);
        intTime = cell2mat(intTimeC');
        expd = fitdist(intTime, 'Exponential');

        % Compute lick counts and fit poisson pdf
        numBins = ceil((minMaxTime(2) - minMaxTime(1)) / binSize);
        binEdges = linspace(minMaxTime(1), minMaxTime(2), numBins + 1);
        counts = cellfun(@(a) histcounts(a, binEdges), timeStampC, 'UniformOutput', false);
        countsMat = cell2mat(counts');
        lambda = mean(countsMat(:));
        poisd = makedist('Poisson', 'lambda', lambda);

        % Simulate using the combined model
        timeStampOut = cell(1, numReps);
        for jj = 1:numReps
            simulatedCounts = random(poisd, numBins, 1);
            timeStampOut{jj} = [];
            for bin = 1:numBins
                binStartTime = binEdges(bin);
                timeStampsInBin = binStartTime + cumsum(random(expd, simulatedCounts(bin), 1));
                timeStampOut{jj} = [timeStampOut{jj}; timeStampsInBin(timeStampsInBin < binEdges(bin+1))];
            end
            intervals = diff([minMaxTime(1); timeStampOut{jj}]);
            intervals = intervals(randperm(length(intervals)));
            timeStampOut{jj} = minMaxTime(1) +cumsum(intervals);
        end
    end

end


function hAx = sigDrezSubplot(plotTitle, sigDrez, subPlotSet, fullColorMap)

hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
hold on;

% Create an empty array to hold the scatter plot handles
h = [];

legendEntries = {};

if isfield(sigDrez{1}, 'hitRate')
    hitRate = cell2mat(cellfun(@(a) a.hitRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(hitRate), hitRate, 100, fullColorMap(1, :), 'filled'); % Append handle
    plot(1:length(hitRate), hitRate, 'Color', fullColorMap(1, :), 'LineStyle', ':');
    legendEntries{end+1} = 'Hit';
end

if isfield(sigDrez{1}, 'missRate')
    missRate = cell2mat(cellfun(@(a) a.missRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(missRate), missRate, 100, fullColorMap(2, :), 'filled'); % Append handle
    plot(1:length(missRate), missRate, 'Color', fullColorMap(2, :), 'LineStyle', ':');
    legendEntries{end+1} = 'Miss';
end

if isfield(sigDrez{1}, 'FaRate')
    FaRate = cell2mat(cellfun(@(a) a.FaRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(FaRate), FaRate, 100, fullColorMap(3, :), 'filled'); % Append handle
    plot(1:length(FaRate), FaRate, 'Color', fullColorMap(3, :), 'LineStyle', ':');
    legendEntries{end+1} = 'FA';
end

if isfield(sigDrez{1}, 'CrRate')
    CrRate = cell2mat(cellfun(@(a) a.CrRate, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(CrRate), CrRate, 100, fullColorMap(4, :), 'filled'); % Append handle
    plot(1:length(CrRate), CrRate, 'Color', fullColorMap(4, :), 'LineStyle', ':');
    legendEntries{end+1} = 'CR';
end

xlim([0.9 length(sigDrez)+0.1]);
ylim([0 1]);

set(gca, 'TickDir', 'out');
xlabel('Blocks');
ylabel('Rate');
title(plotTitle);
hold off;

% Add legend using the scatter plot handles
if ~isempty(legendEntries)
    legend(h, legendEntries, 'Location', 'best');
end

end


function hAx = blockMeanSubPlot(plotTitle, blockMeanC, subPlotSet, fullColorMap, groupNames)
% blockMeanC = {lat.rwdFstLatBlockMean; lat.pnsFstLatBlockMean};
% Use 'cool' colormap and sample it
%cmapIndices = round(linspace(1, size(fullColorMap, 1), length(blockMeanC)));
cmapIndices = round(linspace(1, size(fullColorMap, 1), length(groupNames)));
cmap = fullColorMap(cmapIndices, :);

hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
hold on;
scatterHandles = gobjects(length(blockMeanC), 1); % preallocate a handle array for scatter plots

for i = 1:length(blockMeanC)
    scatterHandles(i) = scatter(1:length(blockMeanC{i}), blockMeanC{i}, 100, cmap(i, :), 'filled');
    plot(1:length(blockMeanC{i}), blockMeanC{i}, 'Color', cmap(i, :), 'LineStyle', ':')
end

% Add the legend using the scatter plot handles and group names
legend(scatterHandles, groupNames, 'Location', 'best'); % 'best' will place the legend at the best location to avoid overlapping with the data.


% Determine the full range of x-values and y-values
allLengths = cellfun(@length, blockMeanC);
allValues = cell2mat(blockMeanC);

minX = 1;
maxX = max(allLengths);
minY = min(allValues(:));
maxY = max(allValues(:));

% Compute 10% margins for x and y axes
xMargin = 0.1 * (maxX - minX);
yMargin = 0.1 * (maxY - minY);

title(plotTitle)
% Update xlim and ylim with the new ranges
xlim([minX - xMargin, maxX + xMargin]);
ylim([minY - yMargin, maxY + yMargin]);

% Set TickDir to 'out', YTick to have six ticks, and YTickLabels with one decimal point
xTicks = 1:length(blockMeanC{i});
yTicks = linspace(minY - yMargin, maxY + yMargin, 6);
yTickLabels = arrayfun(@(x) sprintf('%.1f', x), yTicks, 'UniformOutput', false);

set(gca, 'TickDir', 'out', 'XTick', xTicks, 'YTick', yTicks, 'YTickLabel', yTickLabels) 

hold off;
%print(h, fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101723', 'Figure', 'firstLickLatencyBlocks'), '-dpdf', '-vector');  % '-painters' ensures the output is vector graphics
end


function hAx = sigDrezDprmSubplot(plotTitle, sigDrez, subPlotSet, fullColorMap)

hAx = subplot(subPlotSet(1), subPlotSet(2), subPlotSet(3));
hold on;

% Create an empty array to hold the scatter plot handles
h = [];

legendEntries = {};

if isfield(sigDrez{1}, 'dprime')
    dprime = cell2mat(cellfun(@(a) a.dprime, sigDrez, 'UniformOutput', false));
    h(end+1) = scatter(1:length(dprime), dprime, 100, fullColorMap(6, :), 'filled'); % Append handle
    plot(1:length(dprime), dprime, 'Color', fullColorMap(6, :), 'LineStyle', ':');
    legendEntries{end+1} = 'dprime';
end

xlim([0.9 length(sigDrez)+0.1]);
ylim([0 2.5]);

set(gca, 'TickDir', 'out');
xlabel('Blocks');
ylabel('Rate');
title(plotTitle);
hold off;

% Add legend using the scatter plot handles
if ~isempty(legendEntries)
    legend(h, legendEntries, 'Location', 'best');
end

end




