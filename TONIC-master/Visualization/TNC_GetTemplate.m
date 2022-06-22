function [clstTemp] = TNC_GetTemplate(waveforms,clusterIDs)

uClstIds = unique(clusterIDs);
numClst = numel(uClstIds);

for i=1:numClst
    theseInds = find(clusterIDs==uClstIds(i));
    tmp = waveforms(theseInds,:);
    
    clstTemp.wf(i,:) = mean(tmp,1);
    clstTemp.sd(i,:) = std(tmp,[],1);
    clstTemp.x       = 1:size(waveforms,2);
    
end