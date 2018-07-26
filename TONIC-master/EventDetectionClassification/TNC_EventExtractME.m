%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [spikes] = TNC_EventExtractME(bandPassedData,eventInds,windowPnts)
% FUNCTION DETAILS: This function goes through a single channel of filtered recording data and looks for threshold crossings. A second stage then tests these threshold crossings according to a template matching heuristic to try to classify significant spikes.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS
% spikes.wfs
% spikes.ts
% spikes.tsResolution
% spikes.channel
% spikes.snrThresh
% spikes.winL
% spikes.winR
% spikes.numEvents

debug=0;
    numberOfEvents  = size(eventInds, 2);
    numSamps        = numel(bandPassedData,2);
    numChan         = size(bandPassedData,1);
 
    if numberOfEvents > 0 
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

        if numberOfEvents > 0 
            x = -windowPnts(1,1):windowPnts(1,2);        
            slopes = zeros(numberOfEvents,1);

            % pre-allocate the structure to try and speed this up
            wfs(numberOfEvents).values = zeros(numChan,numel(x));
            if numberOfEvents> 0
                for i=1:numberOfEvents
                    for j=1:numChan
                        mean_ind = eventInds(j,i);       
                        wfs(i).values(j,:) = bandPassedData(j,mean_ind-windowPnts(1,1):mean_ind+windowPnts(1,2));
%                       disp(['i= ' num2str(i) ' j=' num2str(j) ' wfs(i).values(j,:)=' num2str(wfs(i).values(j,:))]);
                    end
                end
            end
        end
    else
        wfs.values = [];
    end
 
% pack the output data structure
    spikes.wfs      = wfs;
    spikes.indsM    = eventInds';
    spikes.inds     = round(mean(eventInds))';  % for Josh
    spikes.x        = x;
    spikes.winL     = windowPnts(1,1);
    spikes.winR     = windowPnts(1,2);
    spikes.numEvs   = numberOfEvents;
    spikes.numChan  = numChan;
    disp(['size(spikes.indsM)=' num2str(size(spikes.indsM)) ' size(spikes.inds)=' num2str(size(spikes.inds))]);
        
% % Plot the detected spikes
%     figure(3); clf;
%     subplot(2,3,1);
%     imagesc(   spikes.wfs   ); hold on;
%     plot([25 25],[0 numberOfEvents],'k--');
%     subplot(2,3,4);
%     imagesc(  corr( spikes.wfs' )  ,[0 1] );
%     subplot(2,3,[2,5]);
%     plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),spikes.wfs,'k');
%     subplot(2,3,[3,6]);
%     plot(1:(windowPnts(1,2)+windowPnts(1,1)+1),spikes.raw,'k');
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
