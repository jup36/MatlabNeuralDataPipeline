function [events] = TNC_EventClean(events);
% FUNCTION DETAILS: a function that provides indices for identified events based upon a generated heuristic function.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS
% requires that EventHeuristic has already been run on the event structure.
% 

    if isfield(events,'heuristic');
        
       heuristicZerosY  = find(events.heuristic.histY==0); % any zeros should be potential threshold positions
       heuristicZerosX  = find(events.heuristic.histX==0); % use this to look for a threshold closest to a zero value
       thresh           = events.heuristic.histX(heuristicZerosY(find(heuristicZerosY>heuristicZerosX,1)));

       % remove events with a heuristic value less than the threshold value,
       % then recreate the histograms
       events.cleanIndices = find(events.heuristic.values > thresh);

    else

        disp('Please run TNC_EventHeuristic before attempting to clean the event data.');

    end


end