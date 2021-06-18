function [ meanW ] = collectAllWeights(weightMat,depths)
    nW = numel(weightMat);
    meanW = nanmean(cell2mat(reshape(weightMat,[1,1,nW])),3); % average across all trial resamples
    meanW(:,4) = depths'./1000; % add the depth in the last column to sort
    meanW(:,5) = 1:size(meanW,1); 
    meanW = sortrows(meanW(~isnan(sum(meanW,2)),:),4); % sort weights by depths
end