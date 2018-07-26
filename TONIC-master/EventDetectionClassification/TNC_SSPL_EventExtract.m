%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [events] = TNC_SSPL_EventExtract(bandPassedData,rawData,eventTimes,windowPnts,dispOn)
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

    numberOfEvents = numel(eventTimes);
    candidateEvents = eventTimes;
 
    numSamples = numel(bandPassedData);
    
    bandPassedData = bandPassedData.*10; % shift to 0.1 uV resolution to allow storing as int_16
    
    waveforms = zeros(numberOfEvents,windowPnts(1,2)+windowPnts(1,1)+1,'int16');
    x = (1:1:size(waveforms,2));
    
%% Check to make sure there is enough room to extract a window
    
    for i=1:numberOfEvents
        
        if candidateEvents(1,i)-windowPnts(1,1) < 1
            
            pad = abs(candidateEvents(1,i)-windowPnts(1,1));
            events.wfs(i).values(1,:) = int16( [zeros(1,pad+1) , bandPassedData(1,1:candidateEvents(1,i)+windowPnts(1,2))] );
            %events.raw(i).values(1,:) = single( [zeros(1,pad+1) , rawData(1,1:candidateEvents(1,i)+windowPnts(1,2))] );     
            
        elseif candidateEvents(1,i)+windowPnts(1,2) > numSamples

            pad = candidateEvents(1,i)+windowPnts(1,2)-numSamples;
            events.wfs(i).values(1,:) = int16( [bandPassedData(1,candidateEvents(1,i)-windowPnts(1,1):numSamples) , zeros(1,pad)] );
            %events.raw(i).values(1,:) = single( [rawData(1,candidateEvents(1,i)-windowPnts(1,1):numSamples) , zeros(1,pad)] );
        
        else
            
            events.wfs(i).values(1,:) = int16( bandPassedData(1,candidateEvents(1,i)-windowPnts(1,1):candidateEvents(1,i)+windowPnts(1,2)) );
            %events.raw(i).values(1,:) = single( rawData(1,candidateEvents(1,i)-windowPnts(1,1):candidateEvents(1,i)+windowPnts(1,2)) );
        
        end
    end
    
% pack the output data structure
%     events.wfs = waveformsFILT;
%     events.raw = waveformsFULL;
    events.inds = candidateEvents;
%     events.ts = candidateEvents./sampleRate;
%     events.tsResolution = sampleRate;
%     events.snrThresh = snrThresh;
    events.winL = windowPnts(1,1);
    events.winR = windowPnts(1,2);
    events.numEvents = numberOfEvents;
    events.resolution = 0.1; %in microvolts;
    
% Plot the detected events
if dispOn
    figure(3); clf;
    subplot(2,3,1);
    imagesc(   events.wfs   ); hold on;
    plot([25 25],[0 numberOfEvents],'k--');
    subplot(2,3,4);
    imagesc(  corr( events.wfs' )  ,[0 1] );
    subplot(2,3,[2,5]);
    plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),events.wfs,'k');
    subplot(2,3,[3,6]);
    plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),events.raw,'k');
    
    drawnow;
end
