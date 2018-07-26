function [meanDataMat] = TNC_GetMeanData(cellArray,errorType);
% FUNCTION DETAILS: Calculate the mean data from a passed cellArray containing matrices of repeated observations (rows) of continuous time series data (each time point is a column). The function returns an object that contains the mean data and a user-definable type of positive and negative error vectors.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% 
% Currently the function supports the following error calculations:
% standard deviation 'std'
% standard error [std / sqrt(N)] 'sem'
% variance 'var'
% full range [maximium and minimum observation at each point] 'range' 
% 
% _________________________________________________________________________
% FREE USE IS PERMITTED AND ENCOURAGED WITH NOTIFICATION OF THE AUTHOR
% 

numChan = size(cellArray);
normCode = 0;

for index = 1:numChan(1,2)
    
    meanDataMat.values(index,:) = mean(cellArray{1,index}(:,:),1);
    
    switch lower(errorType)
        case 'std'
            meanDataMat.posErr(index,:) = std(cellArray{1,index}(:,:),0,1);    
            meanDataMat.negErr(index,:) = std(cellArray{1,index}(:,:),0,1);    
        case 'sem'
            meanDataMat.posErr(index,:) = std(cellArray{1,index}(:,:),0,1);    
            meanDataMat.negErr(index,:) = std(cellArray{1,index}(:,:),0,1);    
            normCode = 1;
        case 'var'
            meanDataMat.posErr(index,:) = var(cellArray{1,index}(:,:),0,1);
            meanDataMat.negErr(index,:) = var(cellArray{1,index}(:,:),0,1);
        case 'range'
            meanDataMat.posErr(index,:) = max(cellArray{1,index}(:,:),[],1);
            meanDataMat.negErr(index,:) = min(cellArray{1,index}(:,:),[],1);
    end
    
end


if normCode
    meanDataMat.posErr = meanDataMat.posErr ./ sqrt(numChan(1,1));
    meanDataMat.negErr = meanDataMat.negErr ./ sqrt(numChan(1,1));
end

meanDataMat.errType = errorType;