function [binSpkMatOut,binEdgesOut] = alignToTaskEvent(evt, spikeTime, binSize, stepSize, oneWindow)           
        binEdgesOut = -oneWindow:binSize:oneWindow; 
        bin1msEdge = cellfun(@(a) max(a-oneWindow,0):a+oneWindow+1, num2cell(evt),'un', 0); 
        bin1msCount = cellfun(@(a) histcounts(spikeTime,a), bin1msEdge, 'un',0); 
        binSpkMatOut = bin1msSpkCountMat(cell2mat(bin1msCount'),binSize,stepSize); 
end