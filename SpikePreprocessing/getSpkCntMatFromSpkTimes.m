function [ spkCntMat ] = getSpkCntMatFromSpkTimes( trBytrSpkTimes, psthParams  )
%This function gets the trial-by-trial spikeCountMat from the trial-by-trial spikeTimes cell 
spkCntMat = cell(size(trBytrSpkTimes,1),1);

if length(psthParams)>1
    psthParams = psthParams{1}; 
end

for i = 1:size(trBytrSpkTimes,1) % increment trial
    spkCntMat{i} = zeros(1,sum(abs(psthParams.psthWin))); % preallocate the spikeCountMat
    
    if isempty(trBytrSpkTimes{i}) % just put 0, when spikeTimes empty
        
    else
        if ~isnan(sum(trBytrSpkTimes{i})) % NaN can be the case for evt out-of-bound case
            spkCntMat{i}(1,trBytrSpkTimes{i}) = 1; % assign 1 to corresponding bins of the spikeTimes
        else 
            spkCntMat{i} = nan(1,sum(abs(psthParams.psthWin))); 
        end
    end
end

end

