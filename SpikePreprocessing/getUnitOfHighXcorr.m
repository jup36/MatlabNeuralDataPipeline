function [ xcorUnitIdx, xcorUnit ] = getUnitOfHighXcorr( unitTimeAcrossTrials, xcorThresholdPer )
%Examine all the neuronal spike-train pairs, and see if any xcorr score surpasses a certain threshold (e.g. 30 %) 
% of the spike count of the fewer spiking unit. 
% At each violation of the threshold, get rid of one unit with a lower spike count 
% of the high-xcorr unit pair.

C = nchoosek(1:size(unitTimeAcrossTrials,1), 2); % all possible pairs

xcor = zeros(size(C,1),1); % cross-correlation of each spike-train pair
xcorUnit = []; % unit with high cross-correlation to be removed 
xcorPair = []; % pair with high cross-correlation to be removed

for pair = 1:size(C,1) % increment unit pairs
    
    if isempty(find(xcorUnit==C(pair,1),1)) && isempty(find(xcorUnit==C(pair,2),1))  % check the cross-correlation unless one of the unit has got caught already
       tempSpikeTrainPair = unitTimeAcrossTrials(C(pair,:),:);      % current unit pair
       xcor(pair) = tempSpikeTrainPair(1,:)*tempSpikeTrainPair(2,:)';         % note that xcorr at zero lag is the dot product of the two spike trains prepared with 1ms bins, i.e., xcorr(tempSpikeTrainPair(1,:),tempSpikeTrainPair(2,:),0)
       [tempMinSumSC, tempMinIdx] = min(sum(tempSpikeTrainPair,2), [], 1); % min of the summed spike counts
       if xcor(pair) > tempMinSumSC*xcorThresholdPer  % if the xcor violates the spike co-occurrence tolerance   
          xcorUnit = [xcorUnit; C(pair,tempMinIdx)];  % collect the xcorUnit
          xcorPair = [xcorPair; C(pair,:)];
       end
    end
    
end

xcorUnitIdx = ones(size(unitTimeAcrossTrials,1),1); % an index for units violating the spike co-occurrence tolerance 
xcorUnitIdx(xcorUnit)=0; % mark the xcorUnit

end

