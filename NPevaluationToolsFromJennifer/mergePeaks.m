function [spikeTimes, spikeAmps, spikeSites] = mergePeaks(allTimes, allAmps, allSites, nSites, hCfg)
    %adapted from JRC, alpha 4.0, commit 4280a02
    %MERGEPEAKS Merge duplicate peak events
    spikeTimes = allTimes';
    spikeAmps = allAmps';
    spikeSites = allSites';

    [spikeTimes, argsort] = sort(spikeTimes);
    spikeAmps = spikeAmps(argsort);
    spikeSites = int32(spikeSites(argsort));
    spikeTimes = int32(spikeTimes);

    [mergedTimes, mergedAmps, mergedSites] = deal(cell(nSites,1));

    try
        parfor iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch ME% don't try to display an error here
               warning('failed to use parallel pool %d: %s', iSite, ME.message);
            end
        end
    catch % parfor failure
        for iSite = 1:nSites
            try
                [mergedTimes{iSite}, mergedAmps{iSite}, mergedSites{iSite}] = ...
                    mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg);
            catch ME
                warning('failed to merge spikes on site %d: %s', iSite, ME.message);
            end
        end
    end

    % merge parfor output and sort
    spikeTimes = neCell2mat(mergedTimes);
    spikeAmps = neCell2mat(mergedAmps);
    spikeSites = neCell2mat(mergedSites);

    [spikeTimes, argsort] = sort(spikeTimes); % sort by time
    spikeAmps = tryGather(spikeAmps(argsort));
    spikeSites = spikeSites(argsort);
end