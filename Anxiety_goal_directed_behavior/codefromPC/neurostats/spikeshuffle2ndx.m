function ndx = spikeshuffle2ndx(ss,L,sortflag,Lstart)
% function ndx = spikeshuffle2ndx(ss,L,sortflag,Lstart)
%
%   using the spikeshuffle structure described in ss, creates a jittered
%    spike train using a partition of size L and a randomly chosen starting
%    point for the partition (or Lstart, if it is provided)
%  
%   ss is a cell array as output by ndx2spikeshuffle
%   there is no error checking on ss
%  
%   L is the width of a jitter window (a positive integer)
%  
%   Lstart is the starting point of the jitter partition
%   Lstart defaults to a (discrete uniform) randomly chosen integer in {1,..,L}
%   
%   if ss = ndx2spikeshuffle(ndx0), then ndx = spikeshuffle2ndx(ss,L)
%    will shuffle the trial labels of spikes in ndx0 within each piece of
%    of the L-partition
%  
%   sortflag = true (default) sorts ndx{t}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(L) ~= 1 || L <= 0 || L ~= round(L)
    error('L must be a positive integer scalar')
end

if nargin < 3 || isempty(sortflag)
    sortflag = true;
else
    sortflag = logical(sortflag);
end
if numel(sortflag) ~= 1
    error('sortflag must be a logical scalar')
end

% pick a random partition
if nargin < 4 || isempty(Lstart)
    Lstart = ceil(rand*L);
end
if numel(Lstart) ~= 1 || Lstart <= 0 || Lstart ~= round(Lstart) || Lstart > L
    error('Lstart must be a positive integer scalar <= L')
end

% data
event = ss{1};
N = size(event,1);

% fixed spikes
minfix = ss{2};
if isinf(minfix)
    Nmax = N;
else
    Nmax = minfix-1;
end

% deal with all fixed spikes
if Nmax < 1
    ndx = event2ndx(event,sortflag,ss{4});
    return
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE JITTER BOUNDARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this partition
bnds = (floor((event(1,1)-Lstart)./L)*L+Lstart:L:floor((event(Nmax,1)-Lstart)./L)*L+Lstart+L).';
% add walls to this partition
bnds = unique([ss{3}; bnds]);
% remove ignorable walls
bnds = bnds(sum(bnds <= event(1,1)):sum(bnds <= event(Nmax,1))+1);
M = numel(bnds)-1;

% there should be at least one event in [bnds(1),bnds(2))
% there should be at least one event in [bnds(M),bnds(M+1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHUFFLE THE TRIALS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over walls
js = 1;
for k = 1:M
    
    st = bnds(k);
    en = bnds(k+1);
    
    % move js (since there is at lest one event in the final case, no need
    %  to test for js out of range)
    while event(js,1) < st
        js = js + 1;
    end
    
    % loop over remaining spikes and exit at the end of this partition
    for j = js:Nmax
        if event(j,1) >= en
            j = j - 1;
            break
        end
    end
    
    % shuffle these indices
    if j >= js
        ej = event(js:j,2);
        ej = ej(randperm(j-js+1));
        event(js:j,2) = ej;
    end

    % move ahead
    js = j+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT TO INDICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndx = event2ndx(event,sortflag,ss{4});
