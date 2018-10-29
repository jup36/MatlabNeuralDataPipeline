function [stillTime] = findJsReadyPt(baselineJsTraj, stillTimeReq, SamplingRate)
%Find the time point when the Js is in position. 
% Algorithm: define the latest timepoint at which the Js sits still for 100
% ms around the time point. Input data represents Js velocity, so find the latest time point around which the absolute Js velocity stays for 0 (or less than 1).  

padLR = stillTimeReq/1000*SamplingRate/2; 
baselineJsVelocity = diff([0 baselineJsTraj]); 

periodicAbsVelSum = nan(length(baselineJsVelocity),1); 

for i = 1:length(baselineJsVelocity)
    if i >= padLR && length(baselineJsVelocity)-i>= padLR && baselineJsTraj(i)<nanmean(baselineJsTraj(1:100))-100
        periodicAbsVelSum(i,1) = sum(abs(baselineJsVelocity(i-padLR+1:i+padLR))); % to find a still point get the periodic Js velocity in a sliding window, also enforce to find a point after the Js movement for placement (just to make sure that the Js has moved from the initial position by a distance of 100)
    end
end

stillPoints = find(periodicAbsVelSum==0); 
[~, minPeriodicAbsVelSum] = min(periodicAbsVelSum); 

if ~isempty(stillPoints)
    stillTime = stillPoints(end); 
elseif min(periodicAbsVelSum) < 20 % if there isn't perfect still point take the lowest velocity point ensuring that it was actually pretty low (lower than 20)
    stillTime = minPeriodicAbsVelSum; 
else
    stillTime = nan; 
end

