function j = mean4spikeshuffle2ndx2jpsth(ss,M,L,Lstart,ndy,My)
% function j = mean4spikeshuffle2ndx2jpsth(ss,M,L,Lstart,ndy,My)
%
% computes the expected value of repeated calls to:
%
%   ndx = spikeshuffle2ndx(ss,L,[],Lstart);
%   j = ndx2jpsth(ndx,M,ndy,My);
%
% use [] to get defaults (for example, use Lstart=[] to average over all
%  partitions)
%
% this is implemented as:
%
%   p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart);
%   b = ndx2bins(ndy,My);
%   j = bins2jpsth(p,b);

j = bins2jpsth(mean4spikeshuffle2ndx2bins(ss,M,L,Lstart),ndx2bins(ndy,My));
