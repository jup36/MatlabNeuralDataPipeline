function rez = binnedLickZTertilesWrapper(var, rez)
%% plot the observed lick counts
lickTimeTertiles = cellfun(@(a) rez.lickTimeC(a), var.tertiles, 'UniformOutput', false);

% observed lick rasters
%rasterPlotCell(lickTimeTertiles, 0.5);

% normalize across all trials
[rez.observedLickZ, var.binEdgesLickZ] = binTimestampsZscore([lickTimeTertiles{:}], var);
var.timeX_binnedLickZ = var.binEdgesLickZ(1:end-1)+ var.binSize/2;
rez.smObservedLickZ = smooth2a(rez.observedLickZ, 0, 3);

% normalize within each tertile
rez.observedLickZEachTertiles  = binTimestampsZscoreEachTertile(lickTimeTertiles, var);
rez.smObservedLickZEachTertiles = smooth2a(rez.observedLickZEachTertiles, 0, 3);

tertiles = divideTrials(size(rez.observedLickZ, 1), 3, 0, 0);

rez.meanAcrossTertileZ = zeros(length(tertiles), size(rez.observedLickZ, 2));
rez.semAcrossTertileZ = zeros(length(tertiles), size(rez.observedLickZ, 2));

rez.meanEachTertileZ = zeros(length(tertiles), size(rez.observedLickZ, 2));
rez.semEachTertileZ = zeros(length(tertiles), size(rez.observedLickZ, 2));

for block = 1:length(tertiles)
    [rez.meanAcrossTertileZ(block, :), ~, rez.semAcrossTertileZ(block, :)] = meanstdsem(rez.smObservedLickZ(tertiles{block}(:), :));
    [rez.meanEachTertileZ(block, :), ~, rez.semEachTertileZ(block, :)] = meanstdsem(rez.smObservedLickZEachTertiles(tertiles{block}(:), :));
end

% plotMeanSem(rez.meanAcrossTertileZ, rez.semAcrossTertileZ, var.timeX_binnedLickZ) % plot the across-trial z-score licks
% plotMeanSem(rez.meanEachTertileZ, rez.semEachTertileZ, var.timeX_binnedLickZ) % plot the across-trial z-score licks


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function zscoreMat = binTimestampsZscoreEachTertile(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numBlocks = length(timeStampC);

        for bl = 1:numBlocks
            zscoreC{bl} = binTimestampsZscore(timeStampC{bl}, var);
        end

        zscoreMat = cell2mat([zscoreC(:)]);

    end

end
