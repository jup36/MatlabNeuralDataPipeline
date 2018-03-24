function binMat = tsp2bin(tsp, t0, t1, Fs)
%% transform tsp format to binary matrix 
% tsp: time of spikes, cell structure, unit-second 
% t0: staring time 
% t1: ending time 
% Fs: sampling frame rate 

numTrial = length(tsp);         %number of trials 
if ~exist('Fs', 'var')
    Fs = 1000; 
end
T = ceil((t1-t0)*Fs);           %number of sample bins 
binMat = false(T, numTrial);    %boolean variable for saving space

for m=1:numTrial
    temp = tsp{m}; 
    %remove spikes outside of interval 
    temp(temp<=t0) = []; 
    temp(temp>t1) = []; 
    
    %change related value into 1
    ind = ceil((temp-t0)*Fs); 
    binMat(ind, m) = true; 
end

binMat = sparse(binMat);  %use sparse matrix to save memory 