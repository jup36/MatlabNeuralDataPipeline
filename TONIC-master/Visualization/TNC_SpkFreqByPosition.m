function [spikeMap] = TNC_SpkFreqByPosition(spiketimes,positions,smthParams);
% FUNCTION DETAILS: 
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI|JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% INPUTS:
% 
% OUTPUTS:

if min(positions.x)<0 | min(positions.y)<0
    positions.x = positions.x - min(positions.x)+1;
    positions.y = positions.y - min(positions.y)+1;
end

% get the extrema of the position data
    extCoords.x(1,1) = min(positions.x);
    extCoords.x(1,2) = max(positions.x);
    extCoords.y(1,1) = min(positions.y);
    extCoords.y(1,2) = max(positions.y);

% create a matrix to hold the full range of position data
    durationMap = zeros( ceil(extCoords.y(1,2)-extCoords.y(1,1)+1) , ceil(extCoords.x(1,2)-extCoords.x(1,1)+1) );
    
% go through entire position data and accumulate samples per location
    for i=1:length(positions.x)
        durationMap(round(positions.y(i)),round(positions.x(i))) = durationMap(round(positions.y(i)),round(positions.x(i)))+1;
    end
    
% convert to time
    durationMap = durationMap./positions.samplePerSec;

% accumulate spikes per position
    spikeMap.count = zeros( ceil(extCoords.y(1,2)-extCoords.y(1,1)+1) , ceil(extCoords.x(1,2)-extCoords.x(1,1)+1) );
    for i=1:length(spiketimes)
        index = find(positions.ts>=spiketimes(i),1);
        spikeMap.spikes.x(i) = round(positions.x(index));
        spikeMap.spikes.y(i) = round(positions.y(index));
        spikeMap.count(spikeMap.spikes.y(i),spikeMap.spikes.x(i)) = spikeMap.count(spikeMap.spikes.y(i),spikeMap.spikes.x(i)) + 1;
    end
    
    spikeMap.spikes.ts = spiketimes';
    
    tmp0 = find(durationMap==0);
    tmpD = find(durationMap);
    
    spikeMap.freq = (spikeMap.count ./ (durationMap+0.001));
    
    if smthParams.sigma ~= 0
        spikeMap.smth.mu = smthParams.mu; % 16 is a reasonable default
        spikeMap.smth.sigma = smthParams.sigma;
        spikeMap.smth.gaussian = TNC_CreateGaussian(spikeMap.smth.mu,spikeMap.smth.sigma,smthParams.mu.*2,1);    
        spikeMap.smth.freq = conv2(spikeMap.freq,spikeMap.smth.gaussian'*spikeMap.smth.gaussian,'same');
    else
        disp('Smoothing is not being used. Set smthParams.sigma to a non-zero value to use smoothing of 2D maps');
        spikeMap.smth.mu = smthParams.mu; % 16 is a reaonable default
        spikeMap.smth.sigma = smthParams.sigma;
        spikeMap.smth.gaussian = 0;    
        spikeMap.smth.freq = spikeMap.freq;
    end
    
% also calculate the isi distribution which is helpful later
    spikeMap.spikes.isi = spikeMap.spikes.ts(2:length(spikeMap.spikes.ts)) - spikeMap.spikes.ts(1:length(spikeMap.spikes.ts)-1);
    spikeMap.spikes.isiH = hist(spikeMap.spikes.isi,0:1:100000);    