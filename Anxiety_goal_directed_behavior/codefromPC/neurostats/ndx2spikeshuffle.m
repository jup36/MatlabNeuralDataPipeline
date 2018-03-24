function ss = ndx2spikeshuffle(ndx,walls,fixndx)
% function ss = ndx2spikeshuffle(ndx,walls,fixndx)
%
% creates a data structure that is useful for spike shuffle jittering
%  
% ndx is a cell array of spike indices (or a single vector of spike indices)
% ss = {event,minfix,walls,T}, where
%      event = [ndx2event(ndx(free spikes),true); ndx2event(ndx(fixed spikes),true)];
%      walls = unique(walls(:));  % using input walls
%      minfix is the index of the first fixed spike in event (and inf if no fixed spikes)
%      T = numel(ndx); % if ndx is cell (and T = 1, otherwise)
%
%   walls is a (possibly empty) list of indices that spikes cannot be jittered
%    through (walls must be shared across trials, unlike jitter)
%   spike indices that land exactly on a wall are assumed to be to the right of
%    the wall
%   a spike is allowed to jitter onto a wall to it's left, but cannot cross it
%   a spike is allowed to jitter 1 bin before a wall to it's right, but cannot
%    land on it
%   if every desired ndx is between 1 and M, use walls = [1,M+1] to constrain jittered
%    spikes to stay between 1 and M
%  
%   fixndx is a cell array of indices into ndx; the default is empty
%   the spike indices indicated by fixndx are not jittered
%   if fixndx is just a vector of indices, then these same indices are used
%     for all elements of ndx
%   +inf can be used to get the last spike index (which is not the largest
%     spike index if s is not sorted)
%   (fixed spikes cannot be walls, unlike jitter)

% deal with no cell input
if ~iscell(ndx)
    ndx = {ndx};
end
T = numel(ndx);

% create event
event = ndx2event(ndx);
[N,tmp2] = size(event);

% create the logical fixndx
% fixndx
if nargin < 3 || isempty(fixndx)
    fixndx = {};
end
fixF = ~isempty(fixndx);

% convert fixndx to cell format
if fixF
    if ~iscell(fixndx)
        fixndxtmp = fixndx;
        fixndx = cell(T,1);
        for t = 1:T
            fixndx{t} = fixndxtmp;
        end
    end
end

% convert inf to last spike time
% convert fixndx to logical
if fixF
    for t = 1:T
        fixndxt = fixndx{t}(:);
        nst = numel(ndx{t});
        if nst == 0
            fixndx{t} = false(0);
        elseif ~islogical(fixndxt)
            fixndxt(isinf(fixndxt)) = nst;
            f = false(nst,1);
            f(fixndxt) = true;
            fixndx{t} = f;
        elseif nst ~= numel(fixndxt)
            error('sizes of ndx and fixndx must match when fixndx is logical')
        end
    end
end

% sort
[event,sortndx] = sortrows(event);

% move the fixed spikes to the end
if fixF
    fixvect = ndx2event(fixndx);
    fixvect = logical(fixvect(sortndx,1));
    if any(fixvect)
        event = [event(~fixvect,:) ; event(fixvect,:)];
        minfix = N - sum(fixvect) + 1;
    else
        minfix = inf;
    end
else
    minfix = inf;
end

% store everything in ss
ss = cell(4,1);
ss{1} = event;
ss{2} = minfix;

if nargin < 2 || isempty(walls)
    ss{3} = [];
elseif iscell(walls)
    error('walls cannot be a cell array')
elseif ~isequal(walls,round(walls))
    error('walls must be integer valued')
else
    ss{3} = unique(walls(:));
end

ss{4} = T;
