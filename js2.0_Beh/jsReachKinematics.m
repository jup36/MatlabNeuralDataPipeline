function [ movKins ] = jsReachKinematics( jsTrajUnfilt, pullThreshold, trialType, mass, sgfiltFramelen )
%This function takes 1-ms sampled joystick position trajectory and detect
% pulles; first, the pull that resulted in threshold crossing
%
% output: ctn; continuous kinematics (velocity, acceleration, force)

jsTrajSteps = jsTrajUnfilt; % Js trajectory in steps 
jsDistCoeff = 2*pi*90/4000; % js movement distance conversion coefficient (4000: the number of total edges(rises/falls) per resolution of the encoder)
pushThresholdmm = 50*jsDistCoeff; % Js push threshold in mm
jsTrajmm = jsTrajSteps.*jsDistCoeff; % joystick trajectory in mm
pullThresholdmm = pullThreshold*jsDistCoeff; % pullthreshold in mm
smJsTraj = smooth(jsTrajmm,20)'; % smoothing 
sgJsTrajmm = sgolayfilt(jsTrajmm,3,sgfiltFramelen); % sg filtered trajectory in mm

jsVelmmps = diff([jsTrajmm(1) jsTrajmm])./(1/1000); % Js velocity (mm/s) computed on the unfiltered jsTraj
smJsVel = smooth(jsVelmmps,20)'; % smoothing
sgJsVel = sgolayfilt(jsVelmmps,3,sgfiltFramelen); % sg filtering

jsAclmmpss = diff([smJsVel(1) smJsVel])./(1/1000); % diff([jsVel(1) jsVel])./(1/1000);     % Js acceleration (mm/s^2) computed on the unfiltered jsTraj
smJsAcl = smooth(jsAclmmpss,20)'; % smoothing
sgJsAcl = sgolayfilt(jsAclmmpss,3,sgfiltFramelen); % sg filtering
stillPts = 10; % 10 ms (to measure the amount of Js movement in a 10 ms sliding window

% get periodic summed absolute velocity within a sliding 10ms window
periodicAbsVelSum = zeros(1,length(jsTrajSteps));
for i = 1:length(smJsVel)
    if i < stillPts
        periodicAbsVelSum(1,i) = sum(abs(smJsVel(i:i+stillPts-1))); % to find a still point get the periodic Js velocity in a sliding window, also enforce to find a point after the Js movement for placement (just to make sure that the Js has moved from the initial position by a distance of 100)
    elseif i >= stillPts
        periodicAbsVelSum(1,i) = sum(abs(smJsVel(i-stillPts+1:i)));
    end
end

movKins.jsTrajSteps = jsTrajSteps; 
movKins.jsDistCoeff = jsDistCoeff; 
movKins.pushThresholdmm = pushThresholdmm;
movKins.jsTrajmm = jsTrajmm;
movKins.smJsTraj = smJsTraj; 
movKins.pullThresholdmm = pullThresholdmm; 
movKins.sgJsTrajmm = sgJsTrajmm;
movKins.periodicAbsVelSum = periodicAbsVelSum; 
movKins.trialType = trialType; 
movKins.jsVelmmps = jsVelmmps; 
movKins.smJsVel = smJsVel; 
movKins.sgJsVel = sgJsVel; 
movKins.jsAclmmpss = jsAclmmpss; 
movKins.smJsAcl = smJsAcl;
movKins.sgJsAcl = sgJsAcl;

switch trialType
    case 'sp' % for a successful pull
        [ movKins.pullStart, movKins.pullStop, movKins.pullMaxVel, movKins.pullMaxVelI, movKins.forceMN, movKins.maxForce, movKins.maxForceI, movKins.netForce ] = detectPull(jsTrajmm, smJsVel, smJsAcl, mass, pullThresholdmm, periodicAbsVelSum); 
%         subplot(2,2,1); plot(jsTrajmm); hold on; plot(movKins.pullStart,jsTrajmm(movKins.pullStart),'*b'); plot(movKins.pullStop,jsTrajmm(movKins.pullStop),'*r'); hold off;
%         subplot(2,2,2); plot(smJsVel); hold on; plot(movKins.pullStart,smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, smJsVel(movKins.pullMaxVelI),'og'); hold off; 
%         subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
%         subplot(2,2,4); plot(periodicAbsVelSum)
        
    case 'ps' % for a push
        [ movKins.pushStart, movKins.pushStop, movKins.pushMaxVel, movKins.pushMaxVelI, movKins.forceMN, movKins.maxForce, movKins.maxForceI, movKins.netForce ] = detectPush(jsTrajmm, smJsVel, smJsAcl, mass, pushThresholdmm, periodicAbsVelSum); 
%         subplot(2,2,1); plot(jsTrajmm); hold on; plot(movKins.pushStart,jsTrajmm(movKins.pushStart),'*b'); plot(movKins.pushStop,jsTrajmm(movKins.pushStop),'*r'); hold off;
%         subplot(2,2,2); plot(smJsVel); hold on; plot(movKins.pushStart,smJsVel(movKins.pushStart),'*b'); plot(movKins.pushStop,smJsVel(movKins.pushStop),'*r'); plot(movKins.pushMaxVelI, smJsVel(movKins.pushMaxVelI),'og'); hold off; 
%         subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
%         subplot(2,2,4); plot(periodicAbsVelSum)
        
    case 'pm' % for a premature pull
        [ movKins.pullStart, movKins.pullStop, movKins.pullMaxVel, movKins.pullMaxVelI, movKins.forceMN, movKins.maxForce, movKins.maxForceI, movKins.netForce ] = detectPull(jsTrajmm, smJsVel, smJsAcl, mass, pullThresholdmm, periodicAbsVelSum); 
%         subplot(2,2,1); plot(jsTrajmm); hold on; plot(movKins.pullStart,jsTrajmm(movKins.pullStart),'*b'); plot(movKins.pullStop,jsTrajmm(movKins.pullStop),'*r'); hold off;
%         subplot(2,2,2); plot(smJsVel); hold on; plot(movKins.pullStart,smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, smJsVel(movKins.pullMaxVelI),'og'); hold off; 
%         subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
%         subplot(2,2,4); plot(periodicAbsVelSum)

    case 'pmpp' % for a premature pull then push
         [ movKins.pull, movKins.push ] = detectPullandPush(jsTrajmm, smJsVel, smJsAcl, mass, pullThresholdmm, pushThresholdmm, periodicAbsVelSum);     
%         subplot(2,2,1); plot(jsTrajmm); hold on; plot(movKins.pullStart,jsTrajmm(movKins.pullStart),'*b'); plot(movKins.pullStop,jsTrajmm(movKins.pullStop),'*r'); plot(movKins.pushStart,jsTrajmm(movKins.pushStart),'*b'); plot(movKins.pushStop,jsTrajmm(movKins.pushStop),'*r'); hold off;
%         subplot(2,2,2); plot(smJsVel); hold on; plot(movKins.pullStart,smJsVel(movKins.pullStart),'*b'); plot(movKins.pullStop,smJsVel(movKins.pullStop),'*r'); plot(movKins.pullMaxVelI, smJsVel(movKins.pullMaxVelI),'og'); plot(movKins.pushStart,smJsVel(movKins.pushStart),'*b'); plot(pushStop,smJsVel(pushStop),'*r'); plot(pushMaxVelI, smJsVel(pushMaxVelI),'og'); hold off; 
%         subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
%         subplot(2,2,4); plot(periodicAbsVelSum)
end

end

