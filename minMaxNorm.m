function [minMaxNormMat] = minMaxNorm(cellByTimeMat)
% minMax normalization of each time series (row)
        for cc = 1:size(cellByTimeMat,1)
            maxV = max(cellByTimeMat(cc,:),[],2); 
            minV = min(cellByTimeMat(cc,:),[],2); 
            minMaxNormMat(cc,:) = (cellByTimeMat(cc,:)-minV)/(maxV-minV);  
        end        
end
