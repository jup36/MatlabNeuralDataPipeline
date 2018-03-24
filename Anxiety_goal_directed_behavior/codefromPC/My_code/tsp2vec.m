function tspAll = tsp2vec(tsp, t0, t1)
%% convert spikes in all trials/cells into one long vector 
numTrial = length(tsp); 
if ~exist('t0', 'var')
    t0 = -inf; 
end

if ~exist('t1', 'var')
    t1 = inf; 
end

total = 0; 
estFR = 30; 
tspAll = zeros(numTrial*(t1-t0)*estFR,1); 
for m=1:numTrial
    temp = tsp{m}; 
    temp(temp<=t0) = []; 
    temp(temp>t1) = []; 
    n = length(temp); 
    tspAll((total+1):(total+n)) = temp; 
    total = total+n; 
end
tspAll = tspAll(1:total); 

