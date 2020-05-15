function [timesOut, ampsOut, sitesOut] = mergeSpikesSite(spikeTimes, spikeAmps, spikeSites, iSite, hCfg)
    %MERGESPIKESSITE Merge spikes in the refractory period
    nLims = int32(abs(hCfg.refracIntSamp));

    % find neighboring spikes
    nearbySites = findNearbySites(hCfg.siteLoc, iSite, hCfg.evtDetectRad); % includes iSite
    spikesBySite = arrayfun(@(jSite) find(spikeSites == jSite), nearbySites, 'UniformOutput', 0);
    timesBySite = arrayfun(@(jSite) spikeTimes(spikesBySite{jSite}), 1:numel(nearbySites), 'UniformOutput', 0);
    ampsBySite = arrayfun(@(jSite) spikeAmps(spikesBySite{jSite}), 1:numel(nearbySites), 'UniformOutput', 0);

    iiSite = (nearbySites == iSite);
    iSpikes = spikesBySite{iiSite};
    iTimes = timesBySite{iiSite};
    iAmps = ampsBySite{iiSite};

    % search over peaks on neighboring sites and in refractory period to
    % see which peaks on this site to keep
    keepMe = true(size(iSpikes));
    for jjSite = 1:numel(spikesBySite)
        jSite = nearbySites(jjSite);

        jSpikes = spikesBySite{jjSite};
        jTimes = timesBySite{jjSite};
        jAmps = ampsBySite{jjSite};

        if iSite == jSite
            delays = [-nLims:-1, 1:nLims]; % skip 0 delay
        else
            delays = -nLims:nLims;
        end

        for iiDelay = 1:numel(delays)
            iDelay = delays(iiDelay);

            [jLocs, iLocs] = ismember(jTimes, iTimes + iDelay);
            if ~any(jLocs)
                continue;
            end

            % jLocs: index into j{Spikes,Times,Amps}
            % iLocs: (1st) corresponding index into i{Spikes,Times,Amps}/keepMe
            jLocs = find(jLocs);
            iLocs = iLocs(jLocs);

            % drop spikes on iSite where spikes on jSite have larger
            % magnitudes
            nearbyLarger = abs(jAmps(jLocs)) > abs(iAmps(iLocs));
            keepMe(iLocs(nearbyLarger)) = 0;

            % flag equal-amplitude nearby spikes
            ampsEqual = (jAmps(jLocs) == iAmps(iLocs));
            if any(ampsEqual)
                if iDelay < 0 % drop spikes on iSite occurring later than those on jSite
                    keepMe(iLocs(ampsEqual)) = 0;
                elseif iDelay == 0 && jSite < iSite % drop spikes on iSite if jSite is of lower index
                    keepMe(iLocs(ampsEqual)) = 0;
                end
            end
        end
    end

    % keep the peak spikes only
    timesOut = iTimes(keepMe);
    ampsOut = iAmps(keepMe);
    sitesOut = repmat(int32(iSite), size(timesOut));
end