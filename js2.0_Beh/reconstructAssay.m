function [trialInfo] = reconstructAssay(trialsCsv, trStartIdx, rwdIdx)
%This function reads in the trialsCsv file, and trStart and rewardIdx, and
% reconstructs the entire assay structure based off of them. This is
% motivated by curtailing of some of the trials.csv files. 

rPos1s = trialsCsv{trialsCsv.set==0,'reach_position_1'}; % reach Position1s in the 1st set 
pullTqs = trialsCsv{trialsCsv.set==0,'pull_torque'}; % pull torques in the 1st set
trialNumb = trialsCsv{trialsCsv.set==0,'trial'}; % trialNumb within each block of the 1st set

% get the reachPosition1 - pullTorque pairs of the 1st set which is (assumed to be) repeated in the next sets 
posTqPairs = []; 
for t = 1:length(find(trialsCsv.set==0))
    if t == 1
        if t ~= length(find(trialsCsv.set==0))
            tempTrialNumb = trialNumb(find((rPos1s==rPos1s(t)&pullTqs==pullTqs(t))==1,1,'last'))+1;
        end
        posTqPairs = [rPos1s(1), pullTqs(1), tempTrialNumb];      
    elseif rPos1s(t)-rPos1s(t-1)~=0 || pullTqs(t)-pullTqs(t-1)~=0
        if t ~= length(find(trialsCsv.set==0))
           tempTrialNumb = trialNumb(find((rPos1s==rPos1s(t)&pullTqs==pullTqs(t))==1,1,'last'))+1;
           if tempTrialNumb==size(trialsCsv,1)
               if ~isempty(posTqPairs(:,3))
                 tempTrialNumb = unique(posTqPairs(:,3)); 
               else 
                 tempTrialNumb = 10; % by default just put 10
               end
           end
        end
        posTqPairs = [posTqPairs; [rPos1s(t), pullTqs(t), tempTrialNumb]]; 
    end
end
clearvars temp* t

% get trial-by-trial reward indices
rwdDetected = histcounts(rwdIdx,trStartIdx);
rwdDetected(1,length(trStartIdx)) = ~isempty(find(rwdIdx>trStartIdx(end),1)); 

% get trial-by-trial reachPos1, pullTorque
cumRwd = 0; 
switchN = 0; 
rcReachPos1 = zeros(length(trStartIdx),1); 
rcPullTq = zeros(length(trStartIdx),1);
for i = 1:length(trStartIdx)
    cumRwd = cumRwd + rwdDetected(i); 
    if i == 1 % 1st trial
        pairId = rem(switchN,size(posTqPairs,1))+1; 
        numbRwdToSwitch = posTqPairs(pairId,3);
        rcReachPos1(i,1) = posTqPairs(pairId,1);
        rcPullTq(i,1) = posTqPairs(pairId,2);
    else 
        if cumRwd == numbRwdToSwitch
           rcReachPos1(i,1) = posTqPairs(pairId,1);
           rcPullTq(i,1) = posTqPairs(pairId,2);
           switchN = switchN+1; % go to the next block
           pairId = rem(switchN,size(posTqPairs,1))+1; % update the pos-torque Id 
           cumRwd = 0; % reset the cumulative # of reward
           numbRwdToSwitch = posTqPairs(pairId,3); 
        else
           rcReachPos1(i,1) = posTqPairs(pairId,1);
           rcPullTq(i,1) = posTqPairs(pairId,2);
        end
    end
end

% get trial-by-trial pullThreshold (assumed to be constant throughout a session)
if length(unique(trialsCsv.pull_threshold))>1
    error('More than one pull_threshold was used in this session?!')
else
    rcPullThreshold(1:length(trStartIdx),1) = unique(trialsCsv.pull_threshold); 
end

trialInfo.reachPos1 = rcReachPos1; 
trialInfo.pullTq = rcPullTq; 
trialInfo.pullThreshold = rcPullThreshold; 

end

