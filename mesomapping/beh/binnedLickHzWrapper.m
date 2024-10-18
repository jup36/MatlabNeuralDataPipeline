function [zrez, zvar] = binnedLickHzWrapper(var, lickTimeC, varargin)
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

%%
% observed lick rasters
%rasterPlotCell(lickTimes, 0.5);

%% count across all trials
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

[zrez.lickHz, zvar.binEdgesLickHz] = binTimestampsHz(lickTimesSort, var);
zrez.smLickHz = smooth2a(zrez.lickHz, 0, 3); % smooth across time

zvar.timeX_binnedLickHz = zvar.binEdgesLickHz(1:end-1)+ var.binSize/2;

%% normalize within each group
zrez.lickGroupHz = binTimestampsGroupHz(lickTimes, var);
zrez.smLickGroupHz = cellfun(@(a) smooth2a(a, 0, 3), zrez.lickGroupHz, 'UniformOutput', false);

%% mean, sem of lickHz, lickGroupZ
zrez.mLickHz = zeros(length(zvar.trGroups), size(zrez.lickHz, 2));
zrez.sLickHz = zeros(length(zvar.trGroups), size(zrez.lickHz, 2));

zrez.mLickGroupHz = zeros(length(zvar.trGroups), size(zrez.lickHz, 2));
zrez.sLickGroupHz = zeros(length(zvar.trGroups), size(zrez.lickHz, 2));

for group = 1:length(zvar.trGroups)
    [zrez.mLickHz(group, :), ~, zrez.sLickHz(group, :)] = meanstdsem(zrez.smLickHz(zvar.trGroups{group}(:), :));
    [zrez.mLickGroupHz(group, :), ~, zrez.sLickGroupHz(group, :)] = meanstdsem(zrez.smLickGroupHz{group});
end


% plotMeanSem(zrez.mLickHz, zrez.sLickHz, zvar.timeX_binnedLickHz) % plot the across-trial z-score licks
% plotMeanSem(zrez.mLickGroupHz, zrez.sLickGroupHz, zvar.timeX_binnedLickHz) % plot the across-trial z-score licks


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function hzC = binTimestampsGroupHz(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numGroups = length(timeStampC);

        for gr = 1:numGroups
            hzC{gr} = binTimestampsHz(timeStampC{gr}, var);
        end

        %zscoreMat = cell2mat([hzC(:)]);

    end

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

%% binTimestampsZscore
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
%% binTimestampsHz
    function [binCountsTrialsHz, binEdges] = binTimestampsHz(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numTrials = length(timeStampC);

        % Define bin edges
        binEdges = var.minMaxTpre(1):var.binSize:var.minMaxTpost(2); % e.g. from -3s to 8s in 100ms (0.1s) increments

        % Initialize storage for binned and z-score normalized lick counts
        binCountsTrials = zeros(length(timeStampC), length(binEdges) - 1);

        for i = 1:numTrials
            % Bin the lick timestamps for the trial
            binCounts = histc(timeStampC{i}, binEdges);
            if ~isempty(binCounts)
                % Exclude the last bin
                binCounts = binCounts(1:end-1);
                binCountsTrials(i, :) = binCounts;
            end
        end

        % convert to Hz
        binCountsTrialsHz = binCountsTrials./var.binSize;

    end
end
