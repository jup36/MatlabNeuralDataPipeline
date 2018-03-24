function event = ndx2event(ndx,sortflag)
% function event = ndx2event(ndx,sortflag=false)
%
% ndx is a cell array of spike indices (or a single vector of spike indices)
%
% event is an Nx2 matrix, where event(:,1) is the list of all the indices
%  in ndx and event(k,2) is the trial of the index for event(k,1)
%
% if sortflag is true (default = false), then the output will be passed
%  through sortrows

T = numel(ndx);

% empty input
if T == 0
    event = zeros(0,2);
    return
end

% single trial input
if ~iscell(ndx)
    event = ones(T,2);
    event(:,1) = ndx(:);
    return
end

% count spikes to avoid dynamic allocation

ndxlen = zeros(T,1);

for t = 1:T
    ndxlen(t) = numel(ndx{t});
end

ndxlen = cumsum(ndxlen);
N = ndxlen(T);

% assign spikes

event = zeros(N,2);
event(1:ndxlen(1),1) = ndx{1}(:);
event(1:ndxlen(1),2) = 1;

for t = 2:T
    
    event(ndxlen(t-1)+1:ndxlen(t),1) = ndx{t}(:);
    event(ndxlen(t-1)+1:ndxlen(t),2) = t;
    
end

% sorting
if nargin >= 2 && ~isempty(sortflag) && sortflag
    event = sortrows(event); 
end
