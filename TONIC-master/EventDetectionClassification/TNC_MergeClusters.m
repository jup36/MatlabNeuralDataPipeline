function [newIDs] = TNC_MergeClusters(oldIDs,mergeList)

numItems = numel(oldIDs);
newIDs = oldIDs;
numMrgs = numel(mergeList);

for i=2:numMrgs
    
    tmpInds = find(oldIDs == mergeList(i));
    newIDs(tmpInds) = mergeList(1);
    
end