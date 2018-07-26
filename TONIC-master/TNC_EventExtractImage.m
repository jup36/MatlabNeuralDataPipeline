function [eventImages] = TNC_EventExtractImage(dataSource,recStamps,chanArray,windowPnts);
% FUNCTION DETAILS: This function goes through a single channel of filtered recording data and looks for threshold crossings. A second stage then tests these threshold crossings according to a template matching heuristic to try to classify significant events.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS
% events.wfs
% events.ts
% events.tsResolution
% events.channel
% events.snrThresh
% events.winL
% events.winR
% events.numEvents


    waveforms = zeros(numberOfEvents,windowPnts(1,2)+windowPnts(1,1)+1);
    x = (1:1:size(waveforms,2))./sampleRate;
    slopes = zeros(numberOfEvents,1);

    for i=1:numberOfEvents
        waveforms(i,:) = bandPassedData(1,candidateEvents(1,i)-windowPnts(1,1):candidateEvents(1,i)+windowPnts(1,2));
        tmpSlope    = polyfit(x,bandPassedData(1,candidateEvents(1,i)-2:candidateEvents(1,i)+2),1);
        slopes(i,:) = abs(tmpSlope(1,1));
    end
    
% pack the output data structure
    events.wfs = waveforms;
    events.slopes = slopes;
    events.inds = candidateEvents;
    events.ts = candidateEvents./sampleRate;
    events.tsResolution = sampleRate;
    % events.channel = 
    events.snrThresh = snrThresh;
    events.winL = windowPnts(1,1);
    events.winR = windowPnts(1,2);
    events.numEvents = numberOfEvents;
    
% Plot the detected events
    figure(3); clf;
    plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),events.wfs,'k');
    drawnow;