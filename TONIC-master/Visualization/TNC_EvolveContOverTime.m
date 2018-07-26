function [waveRep] = TNC_EvolveContOverTime(chDims,data,timeRange,reMap);
% FUNCTION DETAILS: function is designed to create a matrix of all recorded electrode channels per time point. the time series of matrices are stored in the cell array "waveRep"
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

remSize = size(reMap);

timePnts = timeRange(1,1):1:timeRange(1,2);
timeLngth = length(timePnts);
waveRep = cell(1,timeLngth);

for index = 1:timeLngth
    for jndex = 1:chDims(1,1)*chDims(1,2)       
        if (max(remSize)>1)
            dataIndex = reMap(jndex);
        else
            dataIndex = jndex;
        end  
        [i,j] = ind2sub(chDims,dataIndex);
        waveRep{1,index}(i,j) = data(jndex,timePnts(index));
    end
end

