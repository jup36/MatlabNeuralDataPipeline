function [ pullStart, pullStop, pullMaxVel, pullMaxVelI, forceMN, maxForce, maxForceI, netForce ] = detectPull(jsTrajmm, smJsVel, smJsAcl, mass, pullThresholdmm, periodicAbsVelSum)
%This is a helper function to detect a pull (pull) that crosses the pullThreshold.

%stillPts = 10; % 10ms

fstThresCross = find( jsTrajmm <= pullThresholdmm, 1, 'first'); % detect the pull that crossed the pull threshold
% define pullStart - find the last still point before the 1st threshold crossing
lastStillPt = find(periodicAbsVelSum(1:fstThresCross)<50, 1, 'last'); % this corresponds to the last point at which Js moved less than 0.05mm for the retrospective 10ms (50*(1/1000), as the velocity is in the unit of mm/s)
if ~isempty(lastStillPt)
    %if lastStillPt > stillPts
        pullStart = lastStillPt; %-stillPts;
    %else
        %pullStart = lastStillPt;
    %end
else % if there was no still point at all
    [~,pullStart] = min(periodicAbsVelSum(1:fstThresCross));
end

% define pullStop - look for a still point at which Js moved less than 0.05mm for the recent 10 ms
%[~, velVellIdx] = findpeaks(-smJsVel, 'MinPeakProminence',5); % to define pullStop find valleys
[~, velVellIdx] = findpeaks(-smJsVel, 'MinPeakProminence',3); % to define pullStop find valleys
if ~isempty(find(periodicAbsVelSum(fstThresCross:end)<50,1)) % 1st look if there's any still point after the threshold crossing
    pullStop = fstThresCross+find(periodicAbsVelSum(fstThresCross:end)<50,1)-1;
elseif ~isempty(find(periodicAbsVelSum(fstThresCross:end)<100,1)) % 1st look if there's any still point after the threshold crossing
    pullStop = fstThresCross+find(periodicAbsVelSum(fstThresCross:end)<100,1)-1;
elseif ~isempty(find(velVellIdx>fstThresCross,1))
        finalValley = velVellIdx(find(velVellIdx>fstThresCross,1,'last')); % take the final valley (the point at which the pull velocity starts to decrease) as the pullStop
        [~,temppullStop] = min(abs(smJsVel(finalValley:end))); % find the min vel point after the valley
        pullStop = finalValley + temppullStop-1; % pullStop
elseif ~isempty(find(abs(smJsVel(fstThresCross:end))==0,1))
    pullStop = find(abs(smJsVel(fstThresCross:end))==0,1); 
else % in the lack of any still point
     [~,temppullStop] = min(abs(smJsVel(fstThresCross:end))); % find the min vel point after the valley
     pullStop = fstThresCross + temppullStop-1; % pullStop
end

% find the max velocity point
[pullMaxVel,temppullVelMaxI] = min(smJsVel(pullStart:pullStop));  
pullMaxVelI = temppullVelMaxI+pullStart-1; 

% compute force
forceMN = zeros(2,length(smJsAcl)); 
forceMN(1,:) = smJsAcl/1000*mass; % mN  
velZeroPostMaxVI0 = find(smJsVel(pullMaxVelI:end)>=0,1); % 
if ~isempty(velZeroPostMaxVI0)
    velZeroPostMaxVI = pullMaxVelI+velZeroPostMaxVI0-1; 
else
    velZeroPostMaxVI = NaN; 
end

forceMN(2,pullStart:nanmin(velZeroPostMaxVI,pullStop))=1; 
[maxForce, maxForceI0] = nanmin(forceMN(1,pullStart:nanmin(velZeroPostMaxVI,pullStop)));
maxForceI = pullStart+maxForceI0-1; 
netForce(1,:) = cumtrapz(nanmax(forceMN(1,pullStart:nanmin(velZeroPostMaxVI,pullStop)),0)); % net force in the push direction 
netForce(2,:) = cumtrapz(nanmin(forceMN(1,pullStart:nanmin(velZeroPostMaxVI,pullStop)),0)); % net force in the pull direction 

%pullMaxVel  = -max(velVells(velVellIdx>=pullStart & velVellIdx<=pullStop)); % max pull vel
%pullMaxVelI = velVellIdx(velVellIdx>=pullStart & velVellIdx<=pullStop & velVells == -pullMaxVel); % max pull vel time point

% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pullStart,jsTrajmm(pullStart),'*b'); plot(pullStop,jsTrajmm(pullStop),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pullStart,smJsVel(pullStart),'*b'); plot(pullStop,smJsVel(pullStop),'*r'); plot(pullMaxVelI, smJsVel(pullMaxVelI),'og'); hold off; 
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
% subplot(2,2,4); plot(periodicAbsVelSum)

end

