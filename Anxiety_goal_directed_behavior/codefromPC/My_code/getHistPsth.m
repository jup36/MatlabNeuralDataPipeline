function histPsth = getHistPsth(tsp, t0, t1, Fs, smth_wd)
%%use histogram method to compute psth 
%tsp: time of spikes, cell structure 
%t0: starting time 
%t1: ending time 
%Fs: sampling rate 
%smth_wd: smoothing width 

% numTrial = length(tsp); 
%trials to be used 

%number of frequency 
if ~exist('Fs', 'var')
    Fs = 1000; 
end
%smoothing width 
if ~exist('smth_wd', 'var')
    smth_wd = 0.1; 
end

%get psth by smoothing binary matrix 
binMat = tsp2bin(tsp, t0, t1, Fs); 
binMat = full(binMat); 
T = size(binMat, 1); 

histPsth.x = (1:T)/Fs; 
histPsth.y = smooth(mean(binMat,2), ceil(smth_wd*Fs)); 

