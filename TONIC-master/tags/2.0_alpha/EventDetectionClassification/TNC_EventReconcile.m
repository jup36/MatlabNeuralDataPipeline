function [eventListRec] = TNC_EventReconcile(eventList,minEventTime);
% FUNCTION DETAILS: This function takes in events detected on individual channels and seeks to reconcile all events such that at any moment in time only a single event occurs on a given electrode grouping. The default behavior is to take the channel on which the event with the largest amplitude occurs as the 'reference time' for the event.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS
% eventListRec = eventList with offending events removed
% 
% NEVER MIND THIS. JUST GRAB ALL EVENTS ON ALL CHANNELS AND THEN TAKE ALL
% UNIQUE EVENT TRIGGERS. SO THIS FUNCTION JUST HAS TO GRAB UNIQUE TRIGGER
% EVENTS. ON ANY ONE CHANNEL TWO TRIGGER EVENTS ARE SEPARATED BY A MIN
% DISTANCE, BUT NO SUCH FILTER APPLIES ACROSS CHANNELS.


% concatenate detection on all channels and sort

% look for overlapping events

% decision about events being on too many channels simultaneously

% 