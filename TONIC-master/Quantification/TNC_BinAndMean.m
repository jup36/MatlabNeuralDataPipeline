function [binnedData] = TNC_BinAndMean(xValues,yValues,numBins)

binnedData.numBins = numBins;

xProps.minVal = min(xValues);
xProps.maxVal = max(xValues);

binWidth = round( (xProps.maxVal-xProps.minVal) ./ numBins );
binnedData.binWidth = binWidth;

yProps.minVal = min(yValues);
yProps.maxVal = max(yValues);

binnedData.xProps = xProps;
binnedData.yProps = yProps;

for k=1:numBins

    loEnd = xProps.minVal + ( binWidth.* (k-1) );
    hiEnd = xProps.minVal + ( binWidth.* (k  ) );
    validInds = find(xValues<hiEnd & xValues>loEnd);

    currYvals = yValues(validInds);
    binnedData.bins.center(k)   = (hiEnd+loEnd) ./ 2;
    binnedData.bins.centAct(k)  = mean(xValues(validInds));
    binnedData.bins.centSd(k)   = std(xValues(validInds));
    binnedData.bins.avg(k)      = mean(yValues(validInds));
    binnedData.bins.sd(k)       = std(yValues(validInds));
    binnedData.bins.count(k)    = numel(validInds);
    binnedData.bins.sem(k)      = std(yValues(validInds)) ./ numel(validInds);
    
    clear currYvals;
end

