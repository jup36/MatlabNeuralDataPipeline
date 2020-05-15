% from JRclust 4.0.0-alpha.5
function [peaks, amps] = findPeaks(samplesIn, thresh, nneighBelow, fs)
    %FINDPEAKS Find samples which exceed threshold
    %params for calculating amplitudes
    sT = int32(floor((0.25/1000)*fs));  % req time after batch start
    eT = int32(floor((0.4/1000)*fs));    % req time before batch end
    %amplitude is max - min over the span -st to +eT
    
    if isempty(nneighBelow)
        nneighBelow = 1;
    end

    peaks = [];
    amps = [];
    if isempty(samplesIn)
        return;
    end

    %exceedsThresh = jrclust.utils.tryGather(samplesIn < -abs(thresh));
    exceedsThresh = (samplesIn < -abs(thresh));
    peakLocs = find(exceedsThresh);
    if isempty(peakLocs)
        return;
    end

    % allow only peaks far enough from start to get complete waveform
    ind = find(peakLocs < sT + 1);
    if numel(ind) == numel(peakLocs) % these are the only peaks found, none valid
        return;
    else % discard all with times less than sT + 1
        peakLocs(ind) = [];
    end

    % allow only peaks far enough from end to get complete waveform
    ind = find(peakLocs > numel(samplesIn) - eT );
    if numel(ind) == numel(peakLocs) % these are the only peaks found, none valid
        return;
    else % discard all with times greater than eT
        peakLocs(ind) = [];
    end

    peakCenters = samplesIn(peakLocs);
    % take only "peaks" whose sample neighbors are not larger
    peaks = peakLocs(peakCenters <= samplesIn(peakLocs + 1) & peakCenters <= samplesIn(peakLocs - 1));
    if isempty(peaks)
        return;
    end

    % take only peaks who have one or two sample neighbors exceeding threshold
    if nneighBelow == 1
        peaks = peaks(exceedsThresh(peaks - 1) | exceedsThresh(peaks + 1));
    elseif nneighBelow == 2
        peaks = peaks(exceedsThresh(peaks - 1) & exceedsThresh(peaks + 1));
    end
    
    %single point amplitudes
    amps = zeros(1,numel(peaks));
    for i = 1:numel(peaks)
        amps(i) = abs(samplesIn(peaks(i)));
    end
    
% alternate version: max - min over nearest 0.5 msec
%     for i = 1:numel(peaks)
%         amps(i) = max(samplesIn(peaks(i)-sT:peaks(i)+eT)) - min(samplesIn(peaks(i)-sT:peaks(i)+eT));
%     end
end
