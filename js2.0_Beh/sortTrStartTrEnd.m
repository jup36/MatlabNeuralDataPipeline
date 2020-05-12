function [trEndIdx] = sortTrStartTrEnd(trStartIdxNidq, trEndIdxNidq, trialsCsv, trTimeOut)

trEndIdx = nan(size(trStartIdxNidq)); 

% get trial duration from the csv file
cellTrialsCsv = table2cell(trialsCsv); 
trialDurCsv = nan(size(cellTrialsCsv,1),1); 

for t = 1:size(cellTrialsCsv,1)
    tempS = cellTrialsCsv{t,10};  % trial start stored at csv
    if ~strcmp(cellTrialsCsv{t,13},'NULL')
       tempE = cellTrialsCsv{t,13};
    elseif ~strcmp(cellTrialsCsv{t,16},'NULL')
       tempE = cellTrialsCsv{t,16};
    end
    trialDurCsv(t,1) = 24*3600*(datenum(tempE,'yyyy-mm-dd-hh-MM-ss')-datenum(tempS,'yyyy-mm-dd-hh-MM-ss')); 
end

% sort trial start and trial end according to trial start
trEndIdxVal = zeros(length(trStartIdxNidq),1); 
for t = 1:length(trStartIdxNidq)
    nextTrEnd = trEndIdxNidq(find(trEndIdxNidq>trStartIdxNidq(t),1,'first')); 
    if t<length(trStartIdxNidq)
        if nextTrEnd<trStartIdxNidq(1,t+1)
            trEndIdx(t) = nextTrEnd; 
        else
            trEndIdx(t) = trStartIdxNidq(1,t+1)+ min(trTimeOut, trialDurCsv(t))*25000; 
        end
    else
        if ~isempty(nextTrEnd)
            trEndIdx(t) = nextTrEnd; 
        else
            trEndIdx(t) = trStartIdxNidq(1,t+1)+ min(trTimeOut, trialDurCsv(t))*25000; 
        end
    end
end


end