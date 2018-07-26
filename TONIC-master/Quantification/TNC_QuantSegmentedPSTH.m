function [segResponse] = TNC_QuantSegmentedPSTH(data,segmentLength,kernel)
% FUNCTION DETAILS: function implements the presumed method of a class called a 'responseclass'. Idea is that you might derive some classification of responses based upon clustering, trial type, etc. and mean statistics on these subsets can be calculated/extracted for plotting or further analysis.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

%% Determine the number of segments
totalTrials = size(data,1);

numSegs = floor(totalTrials./segmentLength)
segInds = segmentLength:segmentLength:totalTrials
segResponse.kernel = kernel;
segResponse.segmentLength = segmentLength; 

%% Extract segments, sum and smooth (dependent upon convFlag)
for i = 1:length(segInds)
    
    segResponse.psth(i,:) = sum(data(segInds(i)-segmentLength+1:segInds(i),:),1);
    
    if length(kernel)>1
        if i==1
            disp('...running smoothing...');
        end
        segResponse.psthS(i,:) = conv(segResponse.psth(i,:),kernel,'same'); 
    end

%     disp(i);
    
end

    