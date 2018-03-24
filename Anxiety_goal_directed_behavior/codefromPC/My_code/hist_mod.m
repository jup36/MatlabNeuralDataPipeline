function results = hist_mod(lambda2, meanFR, peakPhase, modFactor, freq_osc, N)
if nargin<6
    N = 101; 
end 
if nargin<5
    freq_osc = 20; 
end
if nargin<4
    modFactor = 0.4;
end
if nargin<3
    peakPhase = pi;
end
if nargin<2
    meanFR = 20;
end
if nargin<1
    ihbasis = hist_kernel(20);      %basis vector
    lambda2 = exp(ihbasis*[-10,-2,0,1, -0.4690]');   %self-recovery function
end

%phase modulatory function
oscFun = @(x, modFactor, peakPhase)(1+modFactor*cos(x-peakPhase));

phi = linspace(-pi, pi, N);
phi0 = phi(1);
dphi = diff(phi(1:2));
A = zeros(N-1);

phi_ind = mod(bsxfun(@plus, 1:(N-1), (1:2000)'), N-1);
oscPhase = phi_ind*dphi+phi0; 
t_post = (1:2000)*dphi/(2*pi*freq_osc); 
Fs = round(1/diff(t_post(1:2))); 
len1 = length(lambda2); 
len2 = round(len1*Fs/1000); 

lambda2((end+1):2000) = 1; 
lambda2 = resample(lambda2, len2, len1); 
lambda2 = lambda2(1:2000); 

lambda3 = oscFun(oscPhase, modFactor, peakPhase);
lambda = meanFR*bsxfun(@times, lambda3, lambda2)/Fs; %conditional firing rate at each time point 

p1 = zeros(size(lambda));  %probability of not firing spike until t
p2 = p1; %probability of firing spike at t
p1(1,:) = 1-lambda(1, :);
p2(1,:) = lambda(1, :);
for m=2:size(p1, 1)
    p1(m, :) = p1(m-1, :).*(1-lambda(m, :));
    p2(m, :) = p1(m-1, :).*lambda(m, :);
end

temp = mod(oscPhase+pi, 2*pi)-pi;
phi_ind = ceil((temp-phi0)/dphi);
phi_ind(phi_ind==0) = N-1;
phi_ind(phi_ind>=N) = N-1;
for m=1:(N-1)
    for n=1:(N-1)
        A(n, m) = sum(p2(phi_ind(:, m)==n, m));
    end
end
% A(:,1) = (A(:,1)+A(:, end))/2;
% A(:, end) = A(:, 1);
% A = A(1:N-1, 1:N-1);
A_norm = bsxfun(@times, A, 1./sum(A,1));
[u, d] = eig(A_norm);
[~, ind] = min(abs(diag(d)-1));

results.phi = phi(2:end); 
results.mod = u(:, ind)/mean(u(:, ind));
results.A = A_norm; 














