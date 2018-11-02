function [outputArg1, outputArg2] = jsReachKinematics( jsTrajUnfilt, pullThreshold, trialType, sgfiltFramelen )
%This function takes 1-ms sampled joystick position trajectory and detect
% reaches; first, the reach that resulted in threshold crossing
%
% output: ctn; continuous kinematics (velocity, acceleration, force)

jsDistCoeff = 2*pi*90/4000; % joystick movement distance conversion coefficient (4000: the number of total edges(rises/falls) per resolution of the encoder)
stillPts = 10; % 10 ms (to measure the amount of Js movement in a 10 ms sliding window
pushThresholdmm = 40*jsDistCoeff; % Js push threshold in mm 

jsDistCoeff = 2*pi*90/4000; % joystick movement distance conversion coefficient (4000: the number of total edges(rises/falls) per resolution of the encoder)
jsTrajmm = jsTrajUnfilt.*jsDistCoeff; % joystick trajectory in mm
pullThresholdmm = pullThreshold*jsDistCoeff; % pullthreshold in mm
sgJsTrajmm = sgolayfilt(jsTrajmm,3,sgfiltFramelen); % sg filtered trajectory in mm

jsVel = diff([jsTrajmm(1) jsTrajmm])./(1/1000); % Js velocity (mm/s) computed on the unfiltered jsTraj
smJsVel = smooth(jsVel,20); % smoothing
sgJsVel = sgolayfilt(jsVel,3,sgfiltFramelen); %sg filtering

jsAcl = diff([jsVel(1) jsVel])./(1/1000);     % Js acceleration (mm/s^2) computed on the unfiltered jsTraj
smJsAcl = smooth(jsAcl,20); % smoothing
sgJsAcl = sgolayfilt(jsAcl,3,sgfiltFramelen); % sg filtering

% get periodic summed absolute velocity within a sliding 10ms window
for i = 1:length(smJsVel)
    if i < stillPts
        periodicAbsVelSum(i,1) = sum(abs(smJsVel(i:i+stillPts-1))); % to find a still point get the periodic Js velocity in a sliding window, also enforce to find a point after the Js movement for placement (just to make sure that the Js has moved from the initial position by a distance of 100)
    elseif i >= stillPts
        periodicAbsVelSum(i,1) = sum(abs(smJsVel(i-stillPts+1:i)));
    end
end

switch trialType
    case 'sp' % for a successful pull
        fstThresCross = find( jsTrajmm < pullThresholdmm, 1, 'first'); % detect the reach that crossed the pull threshold
        % define reachStart - find the last still point before the 1st threshold crossing
        reachStart = find(periodicAbsVelSum(1:fstThresCross)<50, 1, 'last'); % this corresponds to the last point at which Js moved less than 0.5mm for the retrospective 10ms
        % define reachStop
        [velVells, velVellIdx] = findpeaks(-smJsVel, 'MinPeakProminence',3); % to define reachStop find valleys
        if ~isempty(find(velVellIdx>fstThresCross,1))
            reachStop = velVellIdx(find(velVellIdx>fstThresCross,1,'last')); % take the final valley as the reachStop 
        else
            reachStop = length(jsVel); % if not found, just take the final point
        end
        % find the max velocity point
        reachMaxVel  = -max(velVells(velVellIdx>=reachStart & velVellIdx<=reachStop)); % max reach vel 
        reachMaxVelI = velVellIdx(velVellIdx>=reachStart & velVellIdx<=reachStop & velVells == -reachMaxVel); % max reach vel time point
        
    case 'ps' % for a push
        fstThresCross = find( jsTrajmm > pushThresholdmm, 1, 'first'); % detect the reach that crossed the push threshold
        % define pushStart - find the last still point before the 1st threshold crossing
        pushStart = find(periodicAbsVelSum(1:fstThresCross)<50, 1, 'last'); % this corresponds to the last point at which Js moved less than 0.5mm for the retrospective 10ms
        % define pushStop
        [velPeaks, velPeakIdx] = findpeaks(smJsVel, 'MinPeakProminence',3); % to define pushStop find valleys
        if ~isempty(find(velPeakIdx>fstThresCross))
            pushStop = velPeakIdx(find(velPeakIdx>fstThresCross,1,'last')); % take the final peak as the reachStop
        else
            pushStop = length(jsVel);
        end
        
        
        
    case 'im' % for a immature pull 
        
        
    case 'impullandpush' % for a immature pull then push
        
        
end




end

