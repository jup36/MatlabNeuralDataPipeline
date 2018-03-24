function Npred=getNpred(p1, p2, deltat)
N = ceil(deltat*1000)-1; 

p1 = reshape(p1, [], 1); 
p2 = reshape(p2, [], 1); 
tmp_p2 = conv(p2, ones(2*N+1, 1)); 
Npred = sum(p1.*tmp_p2((N+1):(end-N))); 