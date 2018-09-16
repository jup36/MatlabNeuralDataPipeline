function [ immobilesum, immobilecount ] = immobiledetect( statebegin, stateshift, motion, time, timeconstant )
%This function calculates the total time spent in immobility in a given
% time period. This ignores time spent in immobility less than half a second. 
% This function also returns total # of immobility (immobilecount). 
% Modified on 03/20/13 to change the immobile cutoff from .5 to 1.0 sec. 
% So, now ignore immobile time less than a second. 

    currentstatemotion = zeros((stateshift - statebegin +1),2);
    currentstatemotion(:,1) = statebegin:1:stateshift;              % column 1: indices within the current state  
    currentstatemotion(:,2) = motion(statebegin:stateshift,1);      % column 2: bv.motion of corresponding indices     
    currentmotionidx = find(currentstatemotion(:,2)==1);            % index # of rows with motions
    currentinterval(:,1) = zeros(size(find(currentstatemotion(:,2)==1),1)+1,1);      % these are the intervals between motion detects, which indicate time in immobility  
    
    if isempty(currentmotionidx) == 0       % This is for the cases where at least one motion is detected
           for i = 1:size(currentmotionidx,1)+1
               if i == 1;
                  currentinterval(i,1) = (time(currentstatemotion(currentmotionidx(i,1),1),1) - time(currentstatemotion(1,1),1))*timeconstant;          % immobile latency for the 1st interval
               elseif i < size(find(currentstatemotion(:,2)==1),1)+1
                  currentinterval(i,1) = (time(currentstatemotion(currentmotionidx(i,1),1),1) - time(currentstatemotion(currentmotionidx(i-1,1),1),1))*timeconstant;
               else
                  currentinterval(i,1) = (time(currentstatemotion(end,1),1) - time(currentstatemotion(currentmotionidx(i-1,1),1),1))*timeconstant;      % immobile latency for the last interval
               end
           end
           immobilecount = length(find(currentinterval(:,1)>=1));               % immobilecount indicates # of immobile times 
           immobilesum = sum(currentinterval(find(currentinterval(:,1)>=1),1)); % immobilesum contains sum of total valid(>=.5 sec) immobile time in the current state   
    else       % this is the case, where no motion is detected at all, given the state latency greater than 1 sec.  
    immobilesum = (time(stateshift,1) - time(statebegin,1))*timeconstant;        % If no motion is detected the entire time period is of immobility 
    immobilecount = 1;      % be careful to use immobile count as the behavioral index, because it is wierd in a way that complete freezing without any motion will have a low score.
    end
clearvars -except immobilesum immobilecount
return

