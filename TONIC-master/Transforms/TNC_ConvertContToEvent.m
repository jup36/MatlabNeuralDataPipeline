function [event] = TNC_ConvertContToEvent(inputData,channelID);
% FUNCTION DETAILS: event pulses are obtained from continuous recordings and stored as small structures
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

event.chID = channelID;
eventCount = 0;
event.class = 'EV';

% assume a well behaved square pulse train
allSupEvents = find(inputData>(0.5.*max(inputData)));

for index = 2:length(allSupEvents)
    if allSupEvents(index) ~= allSupEvents(index-1)+1 
        eventCount = eventCount+1;
        event.indices(eventCount) = allSupEvents(index);
    end
end

event.count = eventCount;