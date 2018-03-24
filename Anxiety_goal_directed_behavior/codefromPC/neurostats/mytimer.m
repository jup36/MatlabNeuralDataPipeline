function est = mytimer(k,n,flag)
% function est = mytimer(k,n,flag)
% 
% returns the estimated time until completion
% if k out of n steps are finished and the last
% tic was called before the first step
%
% writes the result to the screen

t = toc;
est = t*n/k-t;
t = round(100*t)/100;
est = round(100*est)/100;
if ~flag || any(k == round([2.^(0:40) n./(1:10)]))
disp([num2str(k/n*100) '% complete after ' num2str(t) 's with ' num2str(est) 's remaining'])
end

return