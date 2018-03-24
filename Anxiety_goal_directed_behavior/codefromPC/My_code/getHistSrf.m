function histSrf = getHistSrf(tsp, t0, t1, Fs, maxISI)
%% get self-recovery function with histogram method 
% tsp: time of spikes, cell structure
% t0: starting time 
% t1: ending time 
% Fs: sampling rate 
% maxISI: maximum inter spike interval 
if ~exist('Fs', 'var')||isempty(Fs)
    Fs = 1000; 
end

%maximum influenced ISI
if ~exist('maxISI', 'var')||isempty(maxISI)
    maxISI = 100; 
end

%trials to be used  
binMat = full(tsp2bin(tsp, t0, t1, Fs)); 

%get histogram of t-tstar 
srfBin = 1:(maxISI+1); 
tAfterSpk = ttStar(tsp, t0, t1); 
temp = hist(tAfterSpk(binMat), srfBin); 
temp(end) = []; 
srfHist = smooth(temp)/mean(temp); 
srfBin(end) = []; 

histSrf.x = reshape(srfBin, [], 1); 
histSrf.y = reshape(srfHist, [],1); 