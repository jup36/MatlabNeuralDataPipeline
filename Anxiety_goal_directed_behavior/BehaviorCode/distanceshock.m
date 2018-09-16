function [ B2shockdistmat, B3shockdistmat ] = distanceshock( numbrat, numbses, numbtrial, shockcell )
%'distanceshock calculates the distance of each trial from previous and
%   next shock trials. 
%       input: 
%           numbrat: number of rats
%           numbses: number of sessions   
%           numbtrial: number of trials
%           shockcell: the cell containing the shock schedule
%       output: 
%           shockdistmat: 
%               column1: distance from the previous shock trial
%               column2: distance from the next shock trial

B2shockdistmat = nan(numbtrial,2);      

B3shockdistmat = nan(numbtrial,2);      

currentB2shock = shockcell.B2{numbses,numbrat};
currentB3shock = shockcell.B3{numbses,numbrat}; 

% get distance from prev and next shock trials BLOCK2 
for i = 1:numbtrial
    if i <= 5                           % ignore trials 1-5 because at least 5 trials are required to figure out the shock probability of the given block. 
       B2shockdistmat(i,1) = NaN;       % column1: distance from the previous shock     
       B2shockdistmat(i,2) = NaN;       % column2: distance from the next shock 
       
    elseif i <= max(currentB2shock)                                             % From trial 6 to the max shock trial
       tempprevshocktrial = currentB2shock(1,max(find(currentB2shock < i)));
       tempnextshocktrial = currentB2shock(1,min(find(currentB2shock >= i)));  
       B2shockdistmat(i,1) = i - tempprevshocktrial;                            % column1: distance from the previous shock   
       B2shockdistmat(i,2) = tempnextshocktrial - i + 1;                        % column2: distance from the next shock (+1 is to avoid 0) 
       
    else                                                                        % Trials greater than the max shock trial
       tempprevshocktrial = currentB2shock(1,max(find(currentB2shock < i)));
       B2shockdistmat(i,1) = i - tempprevshocktrial;                            % column1: distance from the previous shock           
       B2shockdistmat(i,2) = numbtrial - i + 1;                                 % column2: distance from the trial end (for the trials whose number greater than the last shock trial, substract from the trial end).             
    end
end
clearvars i 

% get distance from prev and next shock trials BLOCK3  
for i = 1:numbtrial
    if i <= 5                           % ignore trials 1-5 because at least 5 trials are required to figure out the shock probability of the given block. 
       B3shockdistmat(i,1) = NaN;       % column1: distance from the previous shock     
       B3shockdistmat(i,2) = NaN;       % column2: distance from the next shock 
       
    elseif i <= max(currentB3shock)                                             % From trial 6 to the max shock trial
       tempprevshocktrial = currentB3shock(1,max(find(currentB3shock < i)));
       tempnextshocktrial = currentB3shock(1,min(find(currentB3shock >= i)));  
       B3shockdistmat(i,1) = i - tempprevshocktrial;                            % column1: distance from the previous shock   
       B3shockdistmat(i,2) = tempnextshocktrial - i + 1;                        % column2: distance from the next shock
       
    else                                                                        % Trials greater than the max shock trial
       tempprevshocktrial = currentB3shock(1,max(find(currentB3shock < i)));
       B3shockdistmat(i,1) = i - tempprevshocktrial;                            % column1: distance from the previous shock           
       B3shockdistmat(i,2) = numbtrial - i + 1;                                 % column2: distance from the trial end (for the trials whose number greater than the last shock trial, substract from the trial end).            
    end
end
clearvars i 
