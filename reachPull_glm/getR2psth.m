function [r2_psth] = getR2psth(evtInMs, timeBin, rateFit, rateDat, binSize, oneWin)

r2 = @(a,b) ones(1,size(a,2))-nansum((a-b).^2)./nansum((a-repmat(nanmean(a,1), size(a,1), 1)).^2); % r-squared = 1-SSres/SStot;
model_psth = func.alignGlmOutToTaskEvent(evtInMs, timeBin, rateFit, binSize, oneWin); % fitted rate aligned to reachStart
spike_psth = func.alignGlmOutToTaskEvent(evtInMs, timeBin, rateDat, binSize, oneWin); % fitted rate aligned to reachStart
r2_psth = max(r2(reshape(spike_psth,[],1),reshape(model_psth,[],1)),0);

end