a = rand(100,1) ;
% put some NaN's
a(randsample(100,20)) = NaN ;
%%interpolate 
x = 1:length(a) ;
a(isnan(a)) = interp1(x(~isnan(a)),a(~isnan(a)),x(isnan(a))) ;
plot(x,a,'.r') 
hold on
plot(x,a,'b')


sideTrjTr100X
sideTrjTr100Y


