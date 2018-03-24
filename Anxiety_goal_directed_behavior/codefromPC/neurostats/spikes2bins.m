function b = spike2bins(x,varargin)
% function b = spikes2bins(s,edges)
% function b = spikes2bins(s,tmin,tstep,tmax)
%
% s is a cell array of spike times
% b(j,k) is the number of spikes in s{k} that are in [edges(j),edges(j+1))
%
% [ndx,M] = spikes2ndx(s,edges);  % or [ndx,M] = spikes2ndx(s,tmin,tstep,tmax);
% b = ndx2bins(ndx,M);

[ndx,M] = spikes2ndx(x,varargin{:});
b = ndx2bins(ndx,M);