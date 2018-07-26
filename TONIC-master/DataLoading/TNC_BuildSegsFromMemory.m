function [PeWin] = TNC_BuildSegsFromMemory(dataMat,eventTimes,window,chArray);
% FUNCTION DETAILS: This function goes through a list of events and retrieves windowed segments from continuous recordings already loaded into memory (generally this is designed to work for smaller file sizes).
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
%
% INPUTS:
% 
% RETURNS: 
% PeWin - a data structure / cell array that contains a matrix corresponding to each electrode channel.
% 

PeWin = cell(1,length(chArray));

% loop through the event times and retrieve the relevant pieces of data
% from the data matrix in memory
for index = 1:length(chArray)

    channel = chArray(index);
    
    for jndex = 1:length(eventTimes)
        PeWin{1,channel}(jndex,:) = dataMat(channel,eventTimes(jndex)-window(1,1):eventTimes(jndex)+window(1,2));
    end
    
end