function [events] = TNC_EventHeuristic(events,centerFreq);
% FUNCTION DETAILS: test the thresholded event data for quality based upon a heuristic that examines the spectral density of the events by projecting events onto a wavelet or sinc function.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS
% adds the 'heuristic' event labels that can be used by subsequent
% processing steps to eliminate some events

tstart = tic;

x = -events.winL:1:events.winR;
template = -sinc(2.*pi.*x.*centerFreq); % use a pulse of frequency centered around 1kHz as a template

for i=1:events.numEvents
    
%     heuristic(i,1) = dot(template',events.wfs(i,:)');
%     heuristic(i,2) = trapz(abs(events.wfs(i,:)));

    heuristic(i,1) = dot(template',events.raw(i,:)') - dot(template',events.wfs(i,:)');
    heuristic(i,2) = trapz( events.raw(i,:) - events.wfs(i,:) );
    
end

events.heuristic.values = heuristic;
events_hist1 = hist(events.heuristic.values(:,1),-2000:50:2000);
events_hist2 = hist(events.heuristic.values(:,2),-10000:200:10000);
% events.heuristic.histY(1,:) = events_hist1;
% events.heuristic.histX = -500:50:2000;

telapsed = toc(tstart);

disp(' ');
disp(sprintf('DETECTION AT SECOND STAGE | Number of Events: %g | Elapsed time (seconds): %g',events.numEvents,telapsed));
disp(' ');

figure(4);
subplot(3,3,1:2)
bar(-10000:200:10000,(events_hist2)); shading flat;
axis([-10000 10000 0 max((events_hist2))]);
axis off;
subplot(3,3,[4,5,7,8])
plot(events.heuristic.values(:,2),events.heuristic.values(:,1),'k.');%,[-5000 10000],[-5000 10000],'r--');
axis([-10000 10000 -2000 2000]);
subplot(3,3,[6,9])
barh(-2000:50:2000,(events_hist1)); shading flat;
axis([0 max((events_hist1)) -5000 5000]);
axis off;