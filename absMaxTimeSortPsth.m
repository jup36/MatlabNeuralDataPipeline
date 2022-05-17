function [absMaxSortMat, sortIdx] = absMaxTimeSortPsth(orgMat, timeBin) 
%This function takes a cell-by-time mat and sort the rows of the matrix based on the peak time of the absolute valued activity  
    [~,mm_maxI] = max(abs(orgMat(:,timeBin)),[],2);
    sort_mm_maxI = sortrows([mm_maxI, (1:length(mm_maxI))'],1); 
    sortIdx = sort_mm_maxI(:,2); 
    absMaxSortMat = orgMat(sortIdx,:); 
end