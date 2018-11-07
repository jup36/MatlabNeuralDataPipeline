function [jsSetTimePt, jsTrajmm, smJsVel, periodicAbsVelSum] = findJsReadyPt2(baselineJsTraj, SamplingRate)
%Find the time point when the Js is in position.
% Algorithm: To detect the time at which Js got in position, this function
% fist detects the stepper movement setting the Js in place, after the
% offset of this Js setting stepper movement find the earliest still point
% where the Js doesn't move 0.05mm for 10 ms.

jsTrajUnfilt = decimate(baselineJsTraj,round(SamplingRate/1000));
jsDistCoeff = 2*pi*90/4000; % joystick movement distance conversion coefficient (4000: the number of total edges(rises/falls) per resolution of the encoder)
stillPts = 10; % 10 ms (to measure the amount of Js movement in a 10 ms sliding window
jsTrajmm = jsTrajUnfilt.*jsDistCoeff; % joystick trajectory in mm
jsVel = diff([jsTrajmm(1) jsTrajmm])./(1/1000); % Js velocity (mm/s) computed on the unfiltered jsTraj
smJsVel = smooth(jsVel,20)'; % smoothing
pullThresholdmm = -3; % 3 mm

% get periodic summed absolute velocity within a sliding 10ms window
periodicAbsVelSum = zeros(1,length(jsTrajUnfilt));
for i = 1:length(smJsVel)
    if i < stillPts
        periodicAbsVelSum(1,i) = sum(abs(smJsVel(i:i+stillPts-1))); % to find a still point get the periodic Js velocity in a sliding window, also enforce to find a point after the Js movement for placement (just to make sure that the Js has moved from the initial position by a distance of 100)
    elseif i >= stillPts
        periodicAbsVelSum(1,i) = sum(abs(smJsVel(i-stillPts+1:i)));
    end
end

% detect the 1st Js setting movement by the apparatus
if ~isempty(find( jsTrajmm <= pullThresholdmm, 1, 'first'))
    fstThresCross = find( jsTrajmm <= pullThresholdmm, 1, 'first'); % detect the jsSetMove that crossed the pull threshold
    % define jsSetMoveStart - find the last still point before the 1st threshold crossing
    lastStillPt = find(periodicAbsVelSum(1:fstThresCross)<100, 1, 'last'); % this corresponds to the last point at which Js moved less than 0.05mm for the retrospective 10ms (100*(1ms/1000ms), as the velocity is in unit of mm/s)
    if ~isempty(lastStillPt)
        if lastStillPt > stillPts
            jsSetMoveStart = lastStillPt-stillPts;
        else
            jsSetMoveStart = lastStillPt;
        end
    else % if there was no still point at all
        [~,jsSetMoveStart] = min(periodicAbsVelSum(1:fstThresCross));
    end
    
    % define jsSetMoveStop
    if ~isempty(find(periodicAbsVelSum(fstThresCross:end)<100,1)) % 1st look if there's any still point after the threshold crossing
        jsSetTime = fstThresCross+find(periodicAbsVelSum(fstThresCross:end)<100,1)-1;
    else
        jsSetTime = nan;
    end
    
else % in rare cases, if there's no Js Set movement (no threshold crossing due to the set movement)
    stillPoints = find(periodicAbsVelSum==0,1,'last');
    if ~isempty(stillPoints)
        jsSetTime = stillPoints;
    elseif min(periodicAbsVelSum) < 100 % if there isn't perfect still point find a point at which there was less than 0.05mm Js movement in the recent 10 ms
        jsSetTime = find(periodicAbsVelSum<100,1,'last');
    else
        jsSetTime = nan;
    end
    
end
jsSetTimePt = jsSetTime*round(SamplingRate/1000); % back to the time point (at 25kHz) from ms

% plot to validate
% subplot(2,2,1); plot(jsTrajmm); hold on; plot(jsSetMoveStart,jsTrajmm(jsSetMoveStart),'*b'); plot(jsSetTime,jsTrajmm(jsSetTime),'*r'); hold off;
% subplot(2,2,2); plot(smJsVel); hold on; plot(jsSetMoveStart,smJsVel(jsSetMoveStart),'*b'); plot(jsSetTime,smJsVel(jsSetTime),'*r'); hold off; 
% subplot(2,2,3); findpeaks(-smJsVel, 'MinPeakProminence',3);
% subplot(2,2,4); plot(periodicAbsVelSum)

end

