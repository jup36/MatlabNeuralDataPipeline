function [finalstatebegin, finalstateshift, finalstatelatency, finalimmobile, finalimmobilecount, finalimmobileITI, finalimmobilecountITI] = statedetect_controlbehavior( indexstate, state, stateITI, entirestate, time, timeconstant, motion )
% This function detects transition of each input entirestate by finding the first
% non-zero row after a given specific input row using the while loop. 
% state is the array of the specific entirestate of interest. 
% Input:
% state: the state #, entirestate: the array of entire state transition, time: the array of entire timestamp, timeconstant: 20/1000 (20 ms), stateITI: the state ITI #          
% Output: 
% immobilepriorITI and immobileinstate have two columns: the 1st column contains absolute time and the 2nd column contains percent time.  

itiarray = find(entirestate(:,1) == stateITI);          % get the point where each state ITI (9,19 or 29) begins

switch isempty(itiarray)
    case 1      % empty trial return NaN for all output variables
        finalstatebegin = NaN; 
        finalstateshift = NaN;
        finalstatelatency = NaN;
        finalimmobile = NaN;
        finalimmobilecount = NaN;
        finalimmobileITI = NaN;
        finalimmobilecountITI = NaN;
    
    case 0      % valid trials proceed to the following script
        indexstatearray = find(entirestate(:,1) == indexstate);
        numbrow = size(indexstatearray,1);                      % use the number of indexstates instead of the states to avoid complications due to reminder cues
        itishift = zeros(size(itiarray,1),1);                   % the last ITI is not registered due to shift to the next block
        
        % this part is to label each element of the statebegin relative to the
        % itiarray because statebegin can occur multiple times within one interval
        % due to the reminder cue.
        % statebegin = zeros(size(find(entirestate(:,1) == state && statebegin(:,1) < length(entirestate)),1),2);
        statebegin(:,1) = find(entirestate(:,1) == state);        % get the point where each state begins
        immobileinstate = zeros(size(statebegin,1),2);            % immobileinstate contains immobile time during the current state, the size must be equal to the statebegin. Column 1: raw, Column 2: percent
        immobilecountinstate = zeros(size(statebegin,1),1);       % immobilecountinstate contains the # of immobile times in the current state.
        itilength = zeros(numbrow-1,1);                           % the last ITI is not registered due to shift to the next block
        immobilepriorITI = zeros(size(statebegin,1),2);           % immobilepriorITI contains immobile time during the previous ITI, the size must be equal to the statebegin. Column 1: raw, Column 2: percent
        immobilecountpriorITI = zeros(size(statebegin,1),1);      % immobilecountpriorITI contains the # of immobile times during the previous ITI.
        prioritibegin = zeros(size(statebegin,1),1);              % prioritibegin contains the points where the prior ITI of the current state initiates
        prioritishift = zeros(size(statebegin,1),1);              % prioritishift contains the points where the prior ITI of the current state finishes (shifts)
        
        if isempty(indexstatearray) || isempty(itiarray) || isempty(statebegin) == 1        % return NaN if there is no trial in this block
            finalstatebegin = NaN;
            finalstateshift = NaN;
            finalstatelatency = NaN;
            finalimmobileITI = NaN;
            finalimmobile = NaN;
        else
            if entirestate(end,1) == state       % this occurs when the session has been terminated within the state, get rid of the last two elements in the statebegin array.  
                statebegin = statebegin(1:end-2,1);
                numbrow = numbrow-1;        % if the session was terminated while the instrumental poke latency, exclude the last trial from analysis
            else
            end
            
            %% this loop is to identify each statebegin with the number of trial, this is required because of the reminder cue.
            for u = 1:size(statebegin,1)
                if u == 1
                    statebegin(u,2) = 1;         % the first one is always the first
                elseif statebegin(u,1) < itiarray(end,1)
                    statebegin(u,2) = find(itiarray(:,1) > statebegin(u,1),1);      % find the first element of the itiarray that are greater than the current ispoke state.
                else
                    statebegin(u,2) = numbrow;       % in this case, the state is the last trial of each block which equals to the numbrow
                end
            end
            
            %% this part is to get the point where each SOI terminates and shifts (stateshift).
            for o = 1:size(statebegin,1)
                if entirestate(statebegin(o,1)+1,1)==0;         % In this case, state shift doesn't occur in the very next row of the bv.state.
                    i = 1;
                    while entirestate(statebegin(o,1)+i,1)==0;
                        stateshift(o,1) = statebegin(o,1)+i;
                        i = i + 1;
                    end
                    stateshift(o,1) = stateshift(o,1) + 1;
                else
                    stateshift(o,1) = statebegin(o,1) + 1;      % In this case, state shift occurs in the very next row of the very next row of the bv.state.
                end
                statelatency(o,1) = (time(stateshift(o,1),1) - time(statebegin(o,1),1))*timeconstant;
                
                %% This part is to calculate the time and counts of immobility (> 1 sec):
                % Ignore immobility less than 1 sec.
                if statelatency(o,1) < 1            % if a state has been shifted in a second, assume there is no meaningful immobility
                    immobileinstate(o,1) = 0;        % 1st column contains immobile time
                    immobileinstate(o,2) = 0;        % 2nd column contains percent relative to the statelatency
                    immobilecountinstate(o,1) = 0;
                else                                % if a state latency is greater than 1 sec.
                    % use the function immobiledetect to get the sum of meaningful immobility in this each state
                    [immobileinstate(o,1), immobilecountinstate(o,1)] = immobiledetect( statebegin(o,1), stateshift(o,1), motion, time, timeconstant );
                    immobileinstate(o,2) = immobileinstate(o,1)/statelatency(o,1)*100;      % 2nd column contains percent relative to the statelatency
                end
                
                %% this part is to get the time spent motionless during the prior ITI
                if statebegin(o,1) < itiarray(1,1)      % this means that the current statebegin is still in the first trial (before the first ITI).
                    immobilepriorITI(o,1) = NaN;         % ITI exists only from the second event, so put NaNs for the first row.
                    immobilepriorITI(o,2) = NaN;
                    itilength(o,1) = NaN;
                    prioritibegin(o,1) = NaN;
                    prioritishift(o,1) = NaN;
                else
                    prioritibegin(o,1) = itiarray(max(find(itiarray(:,1) < statebegin(o,1))),1);        % get the point where the current prior iti starts.
                    if entirestate(prioritibegin(o,1)+1,1)==0;           % this must be always the case because ITIs are always timed.
                        j = 1;
                        while entirestate(prioritibegin(o,1)+j,1)==0
                            prioritishift(o,1) =  prioritibegin(o,1)+j;
                            j = j + 1;
                        end
                        prioritishift(o,1) = prioritishift(o,1) + 1;
                    else                                                 % this is actually redundant because ITIs are always timed, but keep this just in case.
                        prioritishift(o,1) = prioritibegin(o,1) + 1;
                    end
                    itilength(o,1) = (time(prioritishift(o,1),1) - time(prioritibegin(o,1)))*timeconstant;
                    % use the function immobiledetect to get the sum of meaningful immobility in each ITI
                    
                    [immobilepriorITI(o,1),immobilecountpriorITI(o,1)] = immobiledetect( prioritibegin(o,1), prioritishift(o,1), motion, time, timeconstant );
                    immobilepriorITI(o,2) = immobilepriorITI(o,1)/itilength(o,1)*100;      % percent relative to the statelatency
                end
            end
            
            finalstatebegin = zeros(numbrow,1);
            finalstateshift = zeros(numbrow,1);
            finalstatelatency = zeros(numbrow,1);
            finalimmobile = zeros(numbrow,2);           % column 1: total immobile time within the state, column 2: percent immobile time during the state.
            finalimmobilecount = zeros(numbrow,1);
            finalimmobileITI = zeros(numbrow,3);        % column 1: ITI length, column 2: immobile time during the previous ITI, column 3: percent immobile time during the previous ITI.
            finalimmobilecountITI = zeros(numbrow,1);
            label = statebegin(:,2);
            
            for j = 1:numbrow
                finalstatebegin(j,1) = statebegin(find(label(:,1) == j,1),1);
                finalstateshift(j,1) = stateshift(max(find(label(:,1) == j)),1);                        % use max because if there're multiple shifts (due to reminder cues), the last one is where the state finally shifts
                finalstatelatency(j,1) = nansum(statelatency(find(label(:,1) == j),1));                 % add up if there're multiple latencies in any split states.
                finalimmobile(j,1) = nansum(immobileinstate(find(label(:,1) == j),1));                  % add up if there're multiple immobile latencies in any split states.
                finalimmobile(j,2) = finalimmobile(j,1)/finalstatelatency(j,1)*100;                     % second column contains percentage
                finalimmobilecount(j,1) = nansum(immobilecountinstate(find(label(:,1) == j),1)); % add up if there're multiple immobile counts in any split states.
                finalimmobileITI(j,1) = itilength(find(label(:,1) == j,1),1);                             % if there're multiple itilengths in the current state, only choose the first one because they must be all the same.
                finalimmobileITI(j,2) = immobilepriorITI(find(label(:,1) == j,1),1);               % if there're multiple immobilepriorITIs in the current state, only choose the first one because they must be all the same and redundant.
                finalimmobileITI(j,3) = finalimmobileITI(j,2)/finalimmobileITI(j,1)*100;        % second column contains percentage
                finalimmobilecountITI(j,1) = immobilecountpriorITI(find(label(:,1) == j,1),1);     % if there're multiple immobilecountpriorITIs in the current state, only choose the first one because they must be all the same and redundant.
            end
        end
        
end
return