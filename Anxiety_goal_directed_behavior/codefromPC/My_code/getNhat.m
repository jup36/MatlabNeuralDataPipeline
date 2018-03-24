function syncNhat = getNhat(psth1, psth2, deltat)
Fs = 1000; 
T = length(psth1); 
D = round(deltat*Fs); 
N = ceil(T/D)*D; 
if mod(T, D)
    psth1((T+1):N) = 0; 
    psth2((T+1):N) = 0; 
end
y1 = sum(reshape(psth1, D, []),1); 
y2 = sum(reshape(psth2, D, []),1); 

syncNhat = sum(y1.*y2); 