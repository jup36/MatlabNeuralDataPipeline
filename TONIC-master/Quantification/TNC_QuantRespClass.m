function [classResponse] = TNC_QuantRespClass(data,classIds)
% FUNCTION DETAILS: function implements the presumed method of a class called a 'responseclass'. Idea is that you might derive some classification of responses based upon clustering, trial type, etc. and mean statistics on these subsets can be calculated/extracted for plotting or further analysis.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

%% Determine the number of classes
classes     = unique(classIds);
numClasses  = length(classes);


%% For each class extract the data and store within the structure

for i=1:numClasses

    currClass = classes(i);
    
    inds = find(classIds==currClass);
    
    classResponse.class(currClass).psth     = mean(data(inds,:),1);
    classResponse.class(currClass).psthSEM  = std(data(inds,:),0,1) ./ (size(data(inds,:),1)-1);
    classResponse.class(currClass).psthZ    = (classResponse.class(currClass).psth - mean(classResponse.class(currClass).psth)) ./ std(classResponse.class(currClass).psth);
    
end