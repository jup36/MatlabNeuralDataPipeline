function [pullStart, pullStop, pullMaxVel, pullMaxVelI, pushStartCr, pushStopCr, pushMaxVel, pushMaxVelICr] = detectPullandPush(jsTrajmm, smJsVel, pullThresholdmm, pushThresholdmm, periodicAbsVelSum)
%This is a helper function to detect a pull (pull) that crosses the pullThreshold.

%% detect a pull first
[ pullStart, pullStop, pullMaxVel, pullMaxVelI ] = detectPull(jsTrajmm, smJsVel, pullThresholdmm, periodicAbsVelSum); 

% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pullStart,jsTrajmm(pullStart),'*b'); plot(pullStop,jsTrajmm(pullStop),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pullStart,smJsVel(pullStart),'*b'); plot(pullStop,smJsVel(pullStop),'*r'); plot(pullMaxVelI, smJsVel(pullMaxVelI),'og'); hold off; 
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
% subplot(2,2,4); plot(periodicAbsVelSum)

%% then detect a push
smJsVelAfterPull  = smJsVel(pullStop:end); % rest of smJsVel after pull stop from which the push should be detected
jsTrajmmAfterPull = cumsum(smJsVelAfterPull)./1000; % rest of jsTrajmm after pull stop from which the push should be detected

if ~isempty(find(jsTrajmmAfterPull>pushThresholdmm,1))
    [ pushStart, pushStop, pushMaxVel, pushMaxVelI ] = detectPush(jsTrajmmAfterPull, smJsVelAfterPull, pushThresholdmm, periodicAbsVelSum(pullStop:end)); 
else 
    pushStart = nan; 
    pushStop = nan; 
    pushMaxVel = nan; 
    pushMaxVelI = nan;
end

pushStartCr = pullStop+pushStart-1; % time correction
pushStopCr  = pullStop+pushStop-1;  % time correction
pushMaxVelICr  = pullStop+pushMaxVelI-1; % time correction

% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pullStart,jsTrajmm(pullStart),'*b'); plot(pullStop,jsTrajmm(pullStop),'*r'); plot(pushStartCr,jsTrajmm(pushStartCr),'*b'); plot(pushStopCr,jsTrajmm(pushStopCr),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pullStart,smJsVel(pullStart),'*b'); plot(pullStop,smJsVel(pullStop),'*r'); plot(pullMaxVelI, smJsVel(pullMaxVelI),'og'); plot(pushStartCr,smJsVel(pushStartCr),'*b'); plot(pushStopCr,smJsVel(pushStopCr),'*r'); plot(pushMaxVelICr, smJsVel(pushMaxVelICr),'og'); hold off; 
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
% subplot(2,2,4); plot(periodicAbsVelSum)


end

