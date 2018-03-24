function tAfterSpk = ttStar(tsp, t0, t1, n, Fs)
%% calculate the time after previous nth spike 
if ~exist('n', 'var')
    n = 1; 
end
if ~exist('Fs', 'var')
    Fs = 1000; 
end 

numTrial = length(tsp); %number of trials 
T = ceil((t1-t0)*Fs);  %total time T
binMat = full(tsp2bin(tsp, t0, t1)); 

indPrev = cumsum(binMat)-binMat-n+1; %index of previous nth spike 
t = repmat((1:T)', 1, numTrial);   %time 
tStar = 0*t; 
for m=1:numTrial
    tSpk  = find(binMat(:,m));  %time of spikes 
    if isempty(tSpk)
        tStar(:,m) = nan; 
        continue; 
    end
    tmp_ind = indPrev(:,m); 
    tmp_ind(tmp_ind<1) = max(tmp_ind); 
    tStar(:, m) = tSpk(tmp_ind); 
end

tAfterSpk = t-tStar; 
tAfterSpk(tAfterSpk<1) = nan; 