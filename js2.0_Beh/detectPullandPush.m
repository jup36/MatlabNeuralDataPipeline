function [pull, push] = detectPullandPush(jsTrajmm, smJsVel, smJsAcl, mass, pullThresholdmm, pushThresholdmm, periodicAbsVelSum)
%This is a helper function to detect a pull (pull) that crosses the pullThreshold.

%% detect a pull first
[ pull.startI, pull.stopI, pull.maxVel, pull.maxVelI, pull.forceMN, pull.maxForce, pull.maxForceI, pull.netForce ] = detectPull(jsTrajmm, smJsVel, smJsAcl, mass, pullThresholdmm, periodicAbsVelSum);

% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pullStart,jsTrajmm(pullStart),'*b'); plot(pullStop,jsTrajmm(pullStop),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pullStart,smJsVel(pullStart),'*b'); plot(pullStop,smJsVel(pullStop),'*r'); plot(pullMaxVelI, smJsVel(pullMaxVelI),'og'); hold off;
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3);
% subplot(2,2,4); plot(periodicAbsVelSum)

%% then detect a push
smJsVelAfterPull  = smJsVel(pull.stopI:end); % rest of smJsVel after pull stop from which the push should be detected
smJsAclAfterPull  = smJsAcl(pull.stopI:end);
jsTrajmmAfterPull = jsTrajmm(pull.stopI:end); % cumsum(smJsVelAfterPull)./1000; % rest of jsTrajmm after pull stop from which the push should be detected
periodicAbsVelSumAfterPull = periodicAbsVelSum(pull.stopI:end);
normPushThresholdmm = pushThresholdmm+nanmean(jsTrajmmAfterPull(1,1:min(50,length(jsTrajmmAfterPull)))); % normalize the pushThresholdmm
if ~isempty(find(jsTrajmmAfterPull> normPushThresholdmm,1))
    [ push.startI0, push.stopI0, push.maxVel, push.maxVelI0, push.forceMN, push.maxForce, push.maxForceI0, push.netForce ] = detectPushPmpp(jsTrajmmAfterPull, smJsVelAfterPull, smJsAclAfterPull, mass, normPushThresholdmm, periodicAbsVelSumAfterPull);
    push.startI = pull.stopI+push.startI0-1; % time correction
    push.stopI  = pull.stopI+push.stopI0-1;  % time correction
    push.maxVelI  = pull.stopI+push.maxVelI0-1; % time correction
    push.maxForceI = pull.stopI+push.maxForceI0-1; % time correction
else
    push.startI = nan;
    push.stopI = nan;
    push.maxVel = nan;
    push.maxVelI = nan;
    push.maxForceI = nan;
end



% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pullStart,jsTrajmm(pullStart),'*b'); plot(pullStop,jsTrajmm(pullStop),'*r'); plot(pushStartCr,jsTrajmm(pushStartCr),'*b'); plot(pushStopCr,jsTrajmm(pushStopCr),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pullStart,smJsVel(pullStart),'*b'); plot(pullStop,smJsVel(pullStop),'*r'); plot(pullMaxVelI, smJsVel(pullMaxVelI),'og'); plot(pushStartCr,smJsVel(pushStartCr),'*b'); plot(pushStopCr,smJsVel(pushStopCr),'*r'); plot(pushMaxVelICr, smJsVel(pushMaxVelICr),'og'); hold off;
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3);
% subplot(2,2,4); plot(periodicAbsVelSum)


end

