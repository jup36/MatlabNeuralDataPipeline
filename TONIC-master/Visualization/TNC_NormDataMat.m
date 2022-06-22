function [normedData] = TNC_NormDataMat(inputData,refCh,zLogic);
% FUNCTION DETAILS:
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

dims = size(inputData);
 
for index=1:dims(1,1)
    
    normedData(index,:) = inputData(index,:) - inputData(refCh,:);

    if zLogic
        lrMEAN = mean(normedData(index,:),2);
        lrSTD = std(normedData(index,:));
        normedData(index,:) = (normedData(index,:)-lrMEAN)./lrSTD;
    end
    
end

% This data should probably also be converted to Z-scores rather than raw
% values.