function [absMaxSortMat, sortIdx] = absMaxTimeSortPsthNeg(orgMat, timeBin) 
%This function takes a cell-by-time mat and sort the rows of the matrix based on the peak time of the absolute valued activity  
    copyMat = orgMat; 
    copyMat(copyMat>0)=NaN;
    [~,mm_maxI] = max(abs(copyMat(:,timeBin)),[],2);
    sort_mm_maxI = sortrows([mm_maxI, (1:length(mm_maxI))'],1); 
    sortIdx = sort_mm_maxI(:,2); 
    absMaxSortMat = orgMat(sort_mm_maxI(:,2),:); 
    %imagesc(absMaxSortMat)
end
