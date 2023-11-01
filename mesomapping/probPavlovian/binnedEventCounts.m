function [binned_counts, eventsInSec, binEdges] = binnedEventCounts( timeStampsInSec, timeRef, left, right, binSize, stepSize )

edges1ms = (timeRef-left)*1000:1:(timeRef+right)*1000; % edges with 1ms bins
events1msBin = histcounts(timeStampsInSec.*1000, edges1ms); % histogram counts using the 1ms bins
eventsInSec = edges1ms(events1msBin>=1)./1000 - timeRef; 
[binned_counts, binEdges] = bin1msSpkCountMat( events1msBin, binSize.*1000, stepSize.*1000 ); % binning with the specified binSize and stepSize

end