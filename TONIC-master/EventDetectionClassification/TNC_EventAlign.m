function [events] = TNC_EventAlign(spikes,upRatio,dispOn)

% FUNCTION DETAILS: align detected events using interpolation to find the threshold crossing accurately.
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

% if clean==1
%     indSize = size(events.wfs);
%     indices = events.cleanIndices;
%     tmp = size(indices);
%     indSize(1,1) = tmp;
% else
    indSize = size(spikes.wfs);    
% end

x   = 1:1:indSize(1,2);
xx  = 1:(1./upRatio):indSize(1,2);

figure(8); clf;
clear events.Sh_wfs;
count = 1;

splineSpk=zeros(indSize(1,1),numel(xx));


upSize          = numel([-12.*upRatio:24.*upRatio]);
events.Sh_wfs   = zeros(indSize(1,1),upSize);
events.WT_wfs   = zeros(indSize(1,1),upSize);
events.inds     = zeros(indSize(1,1),1);

for i=1:1:indSize(1,1)

    % use spline interpolation to get smoother spikes
    splineSpk(i,:) = spline(x,spikes.wfs(i,:),xx);
        
    % rethreshold...
    tmpMin = min(splineSpk(i,((spikes.winL-1).*upRatio):((spikes.winL+3).*upRatio)));
%     tmpMinCN = min(spikes.wfs(i,((spikes.winL-1).*upRatio):((spikes.winL+3).*upRatio)));
    
    % shift to deal with jitter in peak timing
    newInd  = find(splineSpk(i,((spikes.winL-1).*upRatio):((spikes.winL+3).*upRatio))==tmpMin) + ((spikes.winL-1).*upRatio);
    reSamps = newInd-(12.*upRatio):newInd+(24.*upRatio);

    % calculate the CWT coeffs for the spike
    events.Sh_wfs(count,1:numel(reSamps))   = splineSpk(i,reSamps) - mean(splineSpk(i,1:25));
    events.WT_wfs(count,:)                  = real(cwt(events.Sh_wfs(count,:),numel(reSamps),'sym3'));
    events.inds(count)                      = spikes.inds(i);

    count = count+1;
    
end

events.ShWinL = (12.*upRatio);
events.ShWinR = (24.*upRatio);

events.Sh_x = -events.ShWinL:1:events.ShWinR;

dispSamps = round(rand(1,300) .* size(spikes.wfs,1));

if dispOn == 1
    figure(8); hold on;
    subplot(311);
    plot(1:size(spikes.wfs,2),spikes.wfs(dispSamps,:),'k-');

    subplot(312);
    plot(1:size(events.Sh_wfs,2),events.Sh_wfs(dispSamps,:),'r-');

    subplot(313);
    plot(1:size(events.WT_wfs,2),events.WT_wfs(dispSamps,:),'b-');
end

%% Deprecated code for online display
%     figure(9); clf;
%     plot(x,spikes.wfs(j,:),'k.'); hold on;
%     plot(xx,splineSpk,'r-');
%     plot(xx(newInd),splineSpk(newInd),'b.');
%     
%     newInd
%     ((spikes.winL-1).*upRatio)
%     ((spikes.winR-1).*upRatio)
%     size(splineSpk)


%     tmpMax = max(splineSpk(1,1:((spikes.winL-1).*upRatio)+(5.*upRatio)));
%     if abs(tmpMax) > abs(tmpMin)
%         newInd = find(splineSpk(1,1:((spikes.winL-1).*upRatio)+(5.*upRatio))==tmpMax);
%     else
%     end


%     newInd
%     newInd-((spikes.winL-1).*upRatio)
%     ((spikes.winR+spikes.winL).*upRatio)+newInd
%     spikes.winL+spikes.winR
%     size(1:spikes.winL+spikes.winR-5)

%     size(reSamps)
%     tmp = splineSpk(1,reSamps);
%     size(tmp)
