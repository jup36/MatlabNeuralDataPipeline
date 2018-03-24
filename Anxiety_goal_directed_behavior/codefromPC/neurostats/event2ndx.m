function ndx = event2ndx(event,sortflag,T)
% function ndx = event2ndx(event,sortflag=false,T)
%
% ndx is a cell array of spike indices (or a single vector of spike indices)
%
% event is an Nx2 matrix, where event(:,1) is the list of all the indices
%  in ndx and event(k,2) is the trial of the index for event(k,1)
%
% is sortflag is true (default = false) then ndx will be sorted
%
% if T is provided then ndx will be a cell array of size T, adding empty
%  trials or deleting trials from event as necessary

[N,tmp2] = size(event);

% T check
if nargin < 3 || isempty(T)
    T = [];
elseif numel(T) ~= 1 || T < 0 || T ~= round(T) || isinf(T)
    error('T must be a positive integer scalar')
end

% no trials
if N == 0 || isequal(T,0)
    if isempty(T), T = 0; end
    ndx = cell(T,1);
    return
end

if tmp2 ~= 2, error('event must have 2 columns'), end

% sorting
if nargin < 2 || isempty(sortflag) || ~sortflag
    event = sortrows(event,2);
else
    event = sortrows(event,[2 1]);
end

% checking
if event(1,2) < 1 || ~isequal(event(:,2),round(event(:,2))) || any(isinf(event(:,2)))
    error('event(:,2) must be positive integers')
end

% count the events

Tmax = event(N,2);
if isempty(T), T = Tmax; end

ndxlen = zeros(max(T,Tmax),1);
for k = 1:N
    j = event(k,2);
    ndxlen(j) = ndxlen(j) + 1;
end
ndxlen = cumsum(ndxlen);

% copy into ndx
ndx = cell(T,1);
ndx{1} = event(1:ndxlen(1),1);

for t = 2:T
    ndx{t} = event(ndxlen(t-1)+1:ndxlen(t),1);
end
