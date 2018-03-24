function y = smooth_osc(phi, y0, bw)
%% smooth periodic function 
if ~exist('bw', 'var')
    bw = 0.1; 
end
phi = reshape(phi, 1, []); 
y0 = reshape(y0, 1, []); 
y = y0*0; 
for m=1:length(y0)
    x = mod(phi-phi(m)+pi, 2*pi)-pi; 
    ker = exp(-x.^2/(2*bw^2))/sqrt(2*pi*bw^2); 
    y = y+y0(m)*ker; 
end
% y = y/sum(ker); 
y = y/mean(y); 