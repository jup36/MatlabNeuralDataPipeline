function [N12, sync_tsp1, sync_tsp2]=getNobs(tsp1, tsp2, deltat, t0t1)
if nargin<4
    t0t1 = [-inf, inf]; 
end
if iscell(tsp1)
    num_trial = length(tsp1);
    N12 = 0; 
    sync_tsp1 = cell(num_trial, 1); 
    sync_tsp2 = cell(num_trial, 2); 
    for m=1:num_trial
        [temp, sync_tsp1{m}, sync_tsp2{m}] = getN12(tsp1{m}, tsp2{m}, deltat, t0t1); 
        N12 = N12 + temp; 
    end
else
    [N12, sync_tsp1, sync_tsp2] = getN12(tsp1, tsp2, deltat, t0t1); 
end

function [N12, sync_tsp1, sync_tsp2] = getN12(tsp1, tsp2, deltat, t0t1)
Fs = 1000; 
if nargin>3
    t0 = t0t1(1); 
    t1 = t0t1(2); 
    tsp1(tsp1<=t0) = []; 
    tsp1(tsp1>t1) = []; 
    tsp2(tsp2<=t0) = []; 
    tsp2(tsp2>t1) = []; 
end

deltat = ceil(deltat*Fs)/Fs-1/Fs; 
tsp1 = round(tsp1*Fs)/Fs; 
tsp2 = round(tsp2*Fs)/Fs; 
tsp1 = reshape(tsp1, 1, []); 
tsp2 = reshape(tsp2, [], 1); 

%remove spikes that are too close
min_isi = 0.001; 
ind = [false, diff(tsp1)<min_isi];
while any(ind)
    tsp1(ind) = [];
    ind = [false, diff(tsp1)<=min_isi];
end
ind = [false; diff(tsp2)<min_isi];
while any(ind)
    tsp2(ind) = [];
    ind = [false; diff(tsp2)<=min_isi];
end

difft = bsxfun(@minus, tsp1, tsp2); 
temp = abs(difft)<deltat+0.00001; 
% N12 = sum(sum(temp)); 
[a, b, ~] = find((temp.*cumsum(temp)==1).*(temp.*cumsum(temp, 2)==1)); 
sync_tsp1 = tsp1(b); 
sync_tsp2 = tsp2(a); 
N12 = length(a); 
