function j = mean4jitter2ndx2jpsth(jit,M,L,Lstart,ndy,My)
% function j = mean4jitter2ndx2jpsth(jit,M,L,Lstart,ndy,My)
%
% computes the expected value of repeated calls to:
%
%   ndx = jitter2ndx(jit,L,[],Lstart);
%   j = ndx2jpsth(ndx,M,ndy,My);
%
% use [] to get defaults (for example, use Lstart=[] to average over all
%  partitions)
%
% this is implemented as:
%
%   p = mean4jitter2ndx2bins(jit,M,L,Lstart);
%   b = ndx2bins(ndy,My);
%   j = bins2jpsth(p,b);

j = bins2jpsth(mean4jitter2ndx2bins(jit,M,L,Lstart),ndx2bins(ndy,My));
