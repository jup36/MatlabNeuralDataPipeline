function [mapStructure] = TNC_QuantMap(mapStructure)
% FUNCTION DETAILS: analysis of the autocorrelation (peaks, maximal amp, peak fwhm, peak spacing)
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
% will return various quantifications.
% 
% Quantification of Peak-Valley, xcorr significance, center of mass of
% peaks, fwhm of peaks, distance between multiple peaks, geometry of peaks

%% Operations that work primarily on frequency maps
% generate autocorrelations and countour plots
    mapStructure.smth.auto                         = normxcorr2(mapStructure.smth.freq,mapStructure.smth.freq);
    [mapStructure.smth.contour.freq.values,h]      = imcontour(mapStructure.smth.freq);
    
% use the contour plots to get amplitudes and contour centers
% format of .countour is [values=v;numLines=N] then N pairs of line coords [Nx;Ny] and repeats in that way 
    i=1; count=1;
    while i<length(mapStructure.smth.contour.freq.values)
        
        % get the contour value and number of coordinates
        mapStructure.smth.contour.freq.levels(1,count)  = mapStructure.smth.contour.freq.values(1,i);      
        mapStructure.smth.contour.freq.centers(:,count) = mean(mapStructure.smth.contour.freq.values(:,i+1:i+mapStructure.smth.contour.freq.values(2,i)),2);
        tmpWidth = 0;
        for m= i+1:i+mapStructure.smth.contour.freq.values(2,i)
            tmpWidth(m-i,1) = abs(pdist( [mapStructure.smth.contour.freq.centers(:,count)' ; mapStructure.smth.contour.freq.values(:,m)'] ));
        end
        peakWidths(1,count) = mean(tmpWidth);
        i       = i+mapStructure.smth.contour.freq.values(2,i)+1;
        count   = count+1;

    end

% use the contour plots to get amplitudes and contour centers (NOW FOR THE AUTOCORRELATION FUNCTION)
    figure();
    [mapStructure.smth.contour.auto.values,h]      = imcontour(mapStructure.smth.auto);
    i=1; count=1;
    while i<length(mapStructure.smth.contour.auto.values)
        
        % get the contour value and number of coordinates
        mapStructure.smth.contour.auto.levels(1,count) = mapStructure.smth.contour.auto.values(1,i);      
        mapStructure.smth.contour.auto.centers(:,count) = mean(mapStructure.smth.contour.auto.values(:,i+1:i+mapStructure.smth.contour.auto.values(2,i)),2);
        i       = i+mapStructure.smth.contour.auto.values(2,i)+1;
        count   = count+1;

    end
    
%% Operations that work primarily on frequency maps
% Use the contour data to calculate properties of the peaks: 
% 1. Peak height and width in frequency space
    numFreqLevels = size(mapStructure.smth.contour.freq.levels,2);
    
    findThreshLevel     = find( mapStructure.smth.contour.freq.levels >= (mapStructure.smth.contour.freq.levels(1,numFreqLevels).*0.4) - 0.8, 1 );
    halfmaxSpikeRate    = mapStructure.smth.contour.freq.levels(findThreshLevel);
    tmpThreshInd        = find(mapStructure.smth.contour.freq.levels == halfmaxSpikeRate);
    
    maxLevel        = mapStructure.smth.contour.freq.levels(length(mapStructure.smth.contour.freq.levels));
    mapStructure.smth.contour.freq.supraInds    = tmpThreshInd;
    mapStructure.smth.contour.freq.peakContour  = maxLevel;
    figure(1); clf;
    imagesc(mapStructure.smth.freq);
    hold on;
    plot(mapStructure.smth.contour.freq.centers(1,tmpThreshInd),mapStructure.smth.contour.freq.centers(2,tmpThreshInd),'k.');

% Get stats on the peaks in frequency space
    % heights
        mapStructure.smth.peakHeightsFreq = ones(1,1);
        for o=1:length(tmpThreshInd)          
            mapStructure.smth.peakHeightsFreq(1,o) = mapStructure.smth.freq( round(mapStructure.smth.contour.freq.centers(2,tmpThreshInd(o))), round(mapStructure.smth.contour.freq.centers(1,tmpThreshInd(o))) );
        end
    % widths
        mapStructure.smth.peakWidthsFreq = peakWidths(1,tmpThreshInd).*2;
    % histograms of widths and heights
        mapStructure.smth.peakWidthsFreqHist = hist(mapStructure.smth.peakWidthsFreq,0:1:50);
    % spacings
        distances = pdist(mapStructure.smth.contour.freq.centers(:,mapStructure.smth.contour.freq.supraInds)');
        histTest = hist(distances,0:1:300);
        figure(2); bar(histTest);
        nonZeros = find(histTest);
        hitPeak = nonZeros(1); hitTrough = nonZeros(length(nonZeros)); hitonce=0; neverhit=1;
        for j = 2:length(nonZeros)
            if (nonZeros(j) ~= nonZeros(j-1)+1)
                if hitonce
                    hitTrough   = nonZeros(j-1);
                    break
                end
                hitPeak     = nonZeros(j);
                hitonce = 1; neverhit = 0;
            end
        end
        
        if neverhit | hitPeak>150
           hitPeak = 18;
           hitTrough = nonZeros(length(nonZeros));
        end
        centerSpacings = find(distances>=hitPeak-1 & distances<=hitTrough+1);
        mapStructure.smth.fundSpaceDistFreq(1,1) = mean(distances(centerSpacings));
        mapStructure.smth.fundSpaceDistFreq(1,2) = std(distances(centerSpacings));
        histTest(1,1:4) = 0; % disallow self-spacing over a range of the width of the smoothing kernel (in bins).
        mapStructure.smth.fundSpaceDistFreqHist = histTest;

%% Operations to get peak properties from the contour data to calculate properties of the peaks: 
% 2. Peak spacing (mean and std) and rotation (mean and std).
% 3. Grid completeness (all 6 peaks?)

% for conversion to polar coordinates of suprathreshold contour COMS
% find COMs of contours that are approximately at the 40% of the peak of the normalized correlation
    maxLevel = length(mapStructure.smth.contour.auto.levels);
    topLevel = mapStructure.smth.contour.auto.levels(maxLevel);
    topCenter = mapStructure.smth.contour.auto.centers(:,maxLevel);
    tmpThreshInd = find(mapStructure.smth.contour.auto.levels>=topLevel.*0.4);% & mapStructure.smth.contour.auto.levels<=(topLevel.*0.6)+0.09);
    mapStructure.smth.contour.auto.supraInds = tmpThreshInd;
    figure(3); clf;
    imagesc(mapStructure.smth.auto);
    hold on;
    plot(mapStructure.smth.contour.auto.centers(1,tmpThreshInd),mapStructure.smth.contour.auto.centers(2,tmpThreshInd),'k.');

    [polAngles,polLengths] = cart2pol(mapStructure.smth.contour.auto.centers(2,tmpThreshInd)-topCenter(1,1) , mapStructure.smth.contour.auto.centers(1,tmpThreshInd)-topCenter(2,1));

    distances = pdist(mapStructure.smth.contour.auto.centers(:,tmpThreshInd)');
    mapStructure.smth.contour.auto.pdHist = hist(distances,0:1:300);    
    figure(4); bar(mapStructure.smth.contour.auto.pdHist);

% find the first peak of the radius distribution: find non-zero values and grab the first peak defined as a transition from 0 back to 0
    nonZeros = find(mapStructure.smth.contour.auto.pdHist);
    hitPeak = nonZeros(1); hitTrough = nonZeros(length(nonZeros)); hitonce=0;
    for j = 2:length(nonZeros)
        if (nonZeros(j) ~= nonZeros(j-1)+1)
            if hitonce
                hitTrough   = nonZeros(j-1);
                break
            end
            hitPeak     = nonZeros(j);
            hitonce = 1;
        end
    end

    centerSpacings = find(round(distances)>=hitPeak-1 & round(distances)<=hitTrough+1);
    mapStructure.smth.fundSpaceDistAuto(1,1)    = mean(distances(centerSpacings));
    mapStructure.smth.fundSpaceDistAuto(1,2)    = std(distances(centerSpacings));
    
% get the regularity of the autocorrelation in a circle around the center
% point (a la Mosers)
    mapStructure.smth.contour.auto.polHist = hist(polLengths,0:1:300);
    figure(5); bar(mapStructure.smth.contour.auto.polHist);
    
    nonZeros = find(mapStructure.smth.contour.auto.polHist);
    hitPeak = nonZeros(1); hitTrough = nonZeros(length(nonZeros)); hitonce=0;
    for j = 2:length(nonZeros)
        if (nonZeros(j) ~= nonZeros(j-1)+1)
            if hitonce
                hitTrough   = nonZeros(j-1);
                break
            end
            hitPeak     = nonZeros(j);
            hitonce = 1;
        end
    end
    centerSpacings = find(round(polLengths)>=hitPeak-1 & round(polLengths)<=hitTrough+1);
    figure(5);
    polar(polAngles(centerSpacings),polLengths(centerSpacings),'k.');
    
    mapStructure.smth.fundAnglesAuto            = polAngles(centerSpacings);
    mapStructure.smth.fundRhosAuto              = mean(polLengths(centerSpacings));    
    mapStructure.smth.completenessAuto          = length(centerSpacings)./6;
    
% finally take a circular sampling of the autocorrelation function
    anglesForSampling = -pi:0.1:pi;
    radiusForSampling = mapStructure.smth.fundSpaceDistAuto(1,1) .* ones(1,length(anglesForSampling));
    [pos.x,pos.y] = pol2cart(anglesForSampling,radiusForSampling);
    pos.x
    
    mapStructure.smth.circSampAuto = zeros(1,length(anglesForSampling));
    for n=1:length(pos.x)
        mapStructure.smth.circSampAuto(1,n) = mapStructure.smth.auto(ceil(pos.x(n)+mapStructure.smth.contour.auto.centers(2,maxLevel)),ceil(pos.y(n)+mapStructure.smth.contour.auto.centers(1,maxLevel)));
    end
    mapStructure.smth.circSampAuto(2,:) = anglesForSampling;
    
    figure(21);
    clf; hold on;
    imagesc(mapStructure.smth.auto)
    plot(pos.y+mapStructure.smth.contour.auto.centers(1,maxLevel),pos.x+mapStructure.smth.contour.auto.centers(2,maxLevel),'k.')
    
    mapStructure.smth.entFreq = entropy(mapStructure.smth.freq);
    