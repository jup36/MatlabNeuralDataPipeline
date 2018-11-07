function [stillTime] = findJsReadyPt(baselineJsTraj, stillTimeReqMs, samplingRate)
%Find the time point when the Js is in position. 
% Algorithm: define the latest timepoint at which the Js sits still for 10
% ms around the time point. Input data represents Js velocity, so find the latest time point around which the absolute Js velocity stays for 0 (or less than 1).  
% baselineJsTraj is assumed to be sampled at 25kHz (not downsampled)

jsDistCoeff = 2*pi*90/4000; % joystick movement distance conversion coefficient (4000: the number of total edges(rises/falls) per resolution of the encoder)
baselineJsTrajMm = baselineJsTraj.*jsDistCoeff; 
baselineJsVelocity = smooth(diff([baselineJsTrajMm(1) baselineJsTrajMm]),20)'*(SamplingRate); % mm/sec
periodicAbsVelSum = nan(length(baselineJsVelocity),1); 

for i = 1:length(baselineJsVelocity)
    if i < stillTimeReqMs/1000*samplingRate
        periodicAbsVelSum(i,1) = sum(abs(baselineJsVelocity(i:i+stillTimeReqMs/1000*samplingRate-1))); % to find a still point get the periodic Js velocity in a sliding window
    else
        periodicAbsVelSum(i,1) = sum(abs(baselineJsVelocity(i-stillTimeReqMs/1000*samplingRate+1:i))); % to find a still point get the periodic Js velocity in a sliding window
    end
end

stillPoints = find(periodicAbsVelSum==0,1,'last'); 

if ~isempty(stillPoints)
    stillTime = stillPoints; 
elseif min(periodicAbsVelSum) < 50 % if there isn't perfect still point take the lowest velocity point ensuring that it was actually pretty low (lower than 20)
    stillTime = find(periodicAbsVelSum<20,1,'last'); 
else
    stillTime = nan; 
end

