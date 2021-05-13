function [tagSpkTimeMat50msZmean] = binTrAvgZscore(tagSpkTimeMat1ms,kernel)
%tagSpkTimeMat1ms(:,5000:5000+2) = 0; % to remove light artifact
tagSpkTimeMat1msC = mat2cell(tagSpkTimeMat1ms,ones(size(tagSpkTimeMat1ms,1),1), size(tagSpkTimeMat1ms,2)); % for convolution
tagSpkTimeMat1msCconv = cellfun(@(a) conv(a,kernel,'same'),tagSpkTimeMat1msC,'un',0); % convolution
tagSpkTimeMat1msConv = cell2mat(tagSpkTimeMat1msCconv); % back to mat

tagSpkTimeMat50ms = bin1msSpkCountMat( tagSpkTimeMat1msConv, 50, 50, 'align', 'center' )./(50/1000); % bin
tagSpkTimeMat50msM = nanmean(tagSpkTimeMat50ms); % average across trial
[baseM,baseS] = meanstdsem(tagSpkTimeMat50msM(50:100)'); % baseline stat for normalization
tagSpkTimeMat50msZ = (tagSpkTimeMat50ms-repmat(baseM,size(tagSpkTimeMat50ms,1),size(tagSpkTimeMat50ms,2)))./repmat(baseS,size(tagSpkTimeMat50ms,1),size(tagSpkTimeMat50ms,2)); % z-score
tagSpkTimeMat50msZmean = nanmean(tagSpkTimeMat50msZ ,1); 
end