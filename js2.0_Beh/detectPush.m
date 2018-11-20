function [ pushStart, pushStop, pushMaxVel, pushMaxVelI, forceMN, maxForce, maxForceI, netForce  ] = detectPush(jsTrajmm, smJsVel, smJsAcl, mass, pushThresholdmm, periodicAbsVelSum)
%This is a helper function to detect a push (error) that crosses the pushThreshold.

stillPts = 10; % 10ms

if ~isempty(find( jsTrajmm >= pushThresholdmm, 1, 'first'))
    fstThresCross = find( jsTrajmm >= pushThresholdmm, 1, 'first'); % detect the push that crossed the push threshold
elseif ~isempty(find( jsTrajmm >= pushThresholdmm-1, 1, 'first'))
    fstThresCross = find( jsTrajmm >= pushThresholdmm-1, 1, 'first'); % detect the push that crossed the push threshold
else
    [~,fstThresCross] = max( smJsVel); % detect the push that crossed the push threshold
end

% define pushStart - find the last still point before the 1st threshold crossing
lastStillPt = find(periodicAbsVelSum(1:fstThresCross)<50, 1, 'last'); % this corresponds to the last point at which Js moved less than 0.05mm for the retrospective 10ms (50*(1/1000), as the velocity is in the unit of mm/s)
if ~isempty(lastStillPt)
    if lastStillPt > stillPts
        pushStart = lastStillPt-stillPts;
    else
        pushStart = lastStillPt;
    end
else % if there was no still point at all
    [~,pushStart] = min(periodicAbsVelSum(1:fstThresCross));
end

[~, velPeakIdx] = findpeaks(smJsVel, 'MinPeakProminence',3); % to define pushStop find valleys
% define pushStop - look for a still point at which Js moved less than 0.05mm for the recent 10 ms
if ~isempty(find(periodicAbsVelSum(fstThresCross:end)<50,1)) % 1st look if there's any still point after the threshold crossing
    pushStop = fstThresCross+find(periodicAbsVelSum(fstThresCross:end)<50,1)-1;
elseif ~isempty(find(periodicAbsVelSum(fstThresCross:end)<100,1)) % 1st look if there's any still point after the threshold crossing
    pushStop = fstThresCross+find(periodicAbsVelSum(fstThresCross:end)<100,1)-1;
elseif  ~isempty(find(velPeakIdx>fstThresCross,1))
    finalPeak = velPeakIdx(find(velPeakIdx>fstThresCross,1,'last')); % take the final valley (the point at which the push velocity starts to decrease) as the pushStop
    [~,tempPushStop] = min(abs(smJsVel(finalPeak:end))); % find the min vel point after the valley
    pushStop = finalPeak + tempPushStop-1; % pushStop
    
else
    [~,tempPushStop] = min(abs(smJsVel(fstThresCross:end))); % find the min vel point after the valley
    pushStop = fstThresCross + tempPushStop-1; % pushStop
end

% find the max velocity point
[pushMaxVel,tempPushVelMaxI] = max(smJsVel(pushStart:pushStop));  
pushMaxVelI = tempPushVelMaxI+pushStart-1; 

% compute force
forceMN = zeros(2,length(smJsAcl)); 
forceMN(1,:) = smJsAcl/1000*mass; % mN  
velZeroPostMaxVI0 = find(smJsVel(pushMaxVelI:end)>=0,1); % 
if ~isempty(velZeroPostMaxVI0)
    velZeroPostMaxVI = pushMaxVelI+velZeroPostMaxVI0-1; 
else
    velZeroPostMaxVI = NaN; 
end

forceMN(2,pushStart:nanmin(velZeroPostMaxVI,pushStop))=1; 
[maxForce, maxForceI0] = nanmax(forceMN(1,pushStart:nanmin(velZeroPostMaxVI,pushStop)));
maxForceI = pushStart+maxForceI0-1; 
netForce(1,:) = cumtrapz(nanmax(forceMN(1,pushStart:nanmin(velZeroPostMaxVI,pushStop)),0)); % net force in the push direction 
netForce(2,:) = cumtrapz(nanmin(forceMN(1,pushStart:nanmin(velZeroPostMaxVI,pushStop)),0)); % net force in the pull direction 

%pushMaxVel  = max(velPeaks(velPeakIdx>=pushStart & velPeakIdx<=pushStop)); % max push vel
%pushMaxVelI = velPeakIdx(velPeakIdx>=pushStart & velPeakIdx<=pushStop & velPeaks == pushMaxVel); % max push vel time point

% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pushStart,jsTrajmm(pushStart),'*b'); plot(pushStop,jsTrajmm(pushStop),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pushStart,smJsVel(pushStart),'*b'); plot(pushStop,smJsVel(pushStop),'*r'); plot(pushMaxVelI, smJsVel(pushMaxVelI),'og'); hold off; 
% subplot(2,2,3); findpeaks(smJsVel, 'MinPeakProminence',3); 
% subplot(2,2,4); plot(periodicAbsVelSum)

end

