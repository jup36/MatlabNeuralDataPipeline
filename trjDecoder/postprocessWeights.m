function [wXyzDepthBinMean, wXyzDepthBinMed] = postprocessWeights(wXyzDepth) 
    % bin the x,y,z weights by depth bins
    depthInterval = 0.05; % 50 micron interval
    depthBins = min(wXyzDepth(:,4)):depthInterval:max(wXyzDepth(:,4))+.05; 
    
    [~,~,dI] = histcounts(wXyzDepth(:,4), depthBins); 
    
    wXyzDepthBinMean = nan(length(depthBins),5); 
    wXyzDepthBinMed = nan(length(depthBins),5); 
    
    for d = 1:length(wXyzDepthBinMean)
        wXyzDepthBinMean(d,:) = [nanmean(wXyzDepth(dI==d,:),1), sum(dI==d)]; 
        wXyzDepthBinMed(d,:)  = [nanmedian(wXyzDepth(dI==d,:),1), sum(dI==d)]; 
    end
end