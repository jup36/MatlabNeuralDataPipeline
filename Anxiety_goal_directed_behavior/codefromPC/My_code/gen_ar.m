function y = gen_ar(T, freq)
%% generate data have peak frequency at freq with AR2 process 
%x(t) = phi1*x(t-1)+phi2*x(t-2)+w(t), w(t) is gaussian noise

phi1 = 0.9;
phi2 = -0.8;
N = ceil(freq/0.1667*T/1000);

w = randn(N,1);
y = zeros(N,1);
y(1) = w(1);
y(2) = phi1*w(1)+w(2);
for m=3:N
    y(m) = phi1*w(m-1)+phi2*w(m-2)+w(m);
end

if T*N>2^31
    y = resample(y, ceil(T/N), 1);
    y = y(1:T);
else
    y = resample(y, T, N);
end
y = (y-min(y))/mean(y-min(y));      %min(y) = 0, mean(y)=1