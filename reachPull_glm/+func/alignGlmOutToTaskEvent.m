function [winModelRate] = alignGlmOutToTaskEvent(evt, time_bin, modelRate, binSize, oneWindow)
    windowBins = [floor(-oneWindow/binSize)+1,ceil(oneWindow/binSize)]; 
    binEvt = find(histcounts(evt,time_bin)'); % identify bins per event
    windowBound = arrayfun(@(a) a+windowBins,binEvt,'un',0); 
    for t = 1:length(windowBound)
        if windowBound{t}(1)>=1 && windowBound{t}(2)<=length(modelRate)
            winModelRateC{t} = modelRate(windowBound{t}(1):windowBound{t}(2)); 
        else
            winModelRateC{t} = NaN(sum(abs(windowBins))+1,1); 
        end
    end
    %winModelRateC = cellfun(@(a) modelRate(a(1):a(2)), windowBound, 'un',0); 
    winModelRate = cell2mat(winModelRateC);   
end