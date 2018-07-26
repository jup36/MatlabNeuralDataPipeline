function [events] = TNC_SSPL_EventExtractME(bandPassedData,eventInds,windowPnts)
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

debug=0;
    numberOfEvents  = size(eventInds, 2);
    numSamps        = numel(bandPassedData);
    numChan         = size(bandPassedData,1);
  
    bandPassedData = bandPassedData.*10; % switch to 0.1 uV resolution for minimizing storage overhead
    
    jmin = find(eventInds(:,1) == min(eventInds(:,1))); 
    firstValid = find(eventInds(jmin,:)>windowPnts(1,1),1);
    if firstValid > 1
        len = numel(eventInds(1,:));
        eventInds = eventInds(:,firstValid:len);
    end
    numberOfEvents = numel(eventInds(1,:));   

    jmax = find(eventInds(:,length(eventInds(1,:))) == max(eventInds(:,length(eventInds(1,:)))));
    lastValid       = find(eventInds(jmax,:) < (numSamps-windowPnts(1,2)), 1, 'last');
    if lastValid < numel(eventInds(1,:))
        len = numel(eventInds(1,:));
        eventInds = eventInds(:,1:lastValid);                          
    end
    numberOfEvents  = numel(eventInds(1,:));    
 
    waveforms = zeros(numberOfEvents,windowPnts(1,2)+windowPnts(1,1)+1);
    x = -windowPnts(1,1):windowPnts(1,2);        
    slopes = zeros(numberOfEvents,1);

    % pre-allocate the structure to try and speed this up
    wfs(numberOfEvents).values = zeros(numChan,numel(x),'int16');

    if numberOfEvents> 0
        for i=1:numberOfEvents
            mean_ind = single(eventInds(i));
            wfs(i).values(1:numChan,:) = int16( bandPassedData(1:numChan,mean_ind-windowPnts(1,1):mean_ind+windowPnts(1,2)) );
        end
    end
    
% pack the output data structure
    events.wfs      = wfs;
    events.inds     = eventInds';
    events.x        = x;
    events.winL     = windowPnts(1,1);
    events.winR     = windowPnts(1,2);
    events.numEvs   = numberOfEvents;
    events.numChan  = numChan;
    events.resolution = 0.1;
    
% % Plot the detected events
%     figure(3); clf;
%     subplot(2,3,1);
%     imagesc(   events.wfs   ); hold on;
%     plot([25 25],[0 numberOfEvents],'k--');
%     subplot(2,3,4);
%     imagesc(  corr( events.wfs' )  ,[0 1] );
%     subplot(2,3,[2,5]);
%     plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),events.wfs,'k');
%     subplot(2,3,[3,6]);
%     plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),events.raw,'k');
%     drawnow;

%         if debug==1
%             figure(300); clf;
%             for k=1:numChan
%                 subplot(numChan,1,k);
%                 plot(x, wfs(i).values(k,:),'k');
%                 axis([x(1) x(numel(x)) -1000 300]);
%             end
%             drawnow; pause(0.01);
%         end
