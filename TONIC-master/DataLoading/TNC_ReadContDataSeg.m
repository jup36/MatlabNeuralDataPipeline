function [segments,sampling] = TNC_ReadContDataSeg(tsRange,chList,fname);
% FUNCTION DETAILS: grabs segments of continuous data from large multichannel continuous recording data
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% USE: Blackrock, MAT structures
% REQUIREMENTS: 
% pass tsRange values in seconds
% chList must be a 1xN vector
% tsRange is a (1x2) array
% fname should be a complete path

    tsString = ['e:' num2str(tsRange(1,1)) num2str(tsRange(1,2))];

    for index = 1:length(chList)
        elString = ['e:' num2str(chList(1,index))];
        data = openNSx(fname,'read',elString,tsString,'sec','p:double');
        sampling = data.MetaTags.SampleRes;
        segments(:,index) = data.Data
    end
    
% FOR REFERENCE
% example from the loader: "openNSx('report','read','c:\data\sample.ns5','e:15:30', 't:3:10', 'min', 'p:int16');"
% 
% DEPENDENCIES
% openNSx
