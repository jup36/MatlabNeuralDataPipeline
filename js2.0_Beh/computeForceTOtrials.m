function [ forceMN, netForce ] = computeForceTOtrials(smJsAcl, mass)
%This is a helper function to detect a pull (pull) that crosses the pullThreshold.

% compute force
forceMN = zeros(2,length(smJsAcl)); 
forceMN(1,:) = smJsAcl/1000*mass; % mN  
netForce(1,:) = cumtrapz(nanmax(forceMN(1,:),0)); % net force in the push direction 
netForce(2,:) = cumtrapz(nanmin(forceMN(1,:),0)); % net force in the pull direction 

%pullMaxVel  = -max(velVells(velVellIdx>=pullStart & velVellIdx<=pullStop)); % max pull vel
%pullMaxVelI = velVellIdx(velVellIdx>=pullStart & velVellIdx<=pullStop & velVells == -pullMaxVel); % max pull vel time point

% subplot(2,2,1); plot(jsTrajmm); hold on; plot(pullStart,jsTrajmm(pullStart),'*b'); plot(pullStop,jsTrajmm(pullStop),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(pullStart,smJsVel(pullStart),'*b'); plot(pullStop,smJsVel(pullStop),'*r'); plot(pullMaxVelI, smJsVel(pullMaxVelI),'og'); hold off; 
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3); 
% subplot(2,2,4); plot(periodicAbsVelSum)

end