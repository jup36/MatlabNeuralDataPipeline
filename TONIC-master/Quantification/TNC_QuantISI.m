function [isi] = TNC_QuantISI(spiketimes)
% FUNCTION DETAILS: Calculate the properties of the interspike interval distribution. Includes: calculation of instantaneous ISI, creation of histograms with linear and log spacing, projection onto a classifier space, calculation of 1st moment properties
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

sizeOfTimes = size(spiketimes);

% find non-singleton dimension
if sizeOfTimes(1,1) > sizeOfTimes(1,2)
    lengthOfTimes = sizeOfTimes(1,1);
else
    lengthOfTimes = sizeOfTimes(1,2);    
end

% extract the interspike intervals for each event
isi.instant = spiketimes(2:lengthOfTimes) - spiketimes(1:(lengthOfTimes-1));

% store corresponding timestamps for ease of use later
isi.ts = spiketimes(2:lengthOfTimes);

% create histogram with linear spaced array (in ms)
isi.hist.linTimes = 1:1:2e3;
isi.hist.linCount = hist(isi.instant,isi.hist.linTimes);
isi.hist.linCount = isi.hist.linCount ./ trapz(isi.hist.linCount);

% create histogram with log spaced array (in ms)
logSpace          = -0.6:0.1:6;
isi.hist.logTimes = 10.^logSpace;
isi.hist.logCount = hist(isi.instant,isi.hist.logTimes);
isi.hist.logCount = isi.hist.logCount ./ trapz(isi.hist.logCount);

% calculate stats of first moment of distribution
isi.stats.mean  = mean(isi.instant);
isi.stats.median= median(isi.instant);
isi.stats.var   = var(isi.instant);
isi.stats.std   = std(isi.instant);
isi.stats.fano  = isi.stats.var ./ isi.stats.mean;
