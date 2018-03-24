function s = jitter2ndx(jit,L,sortflag,Lstart)
% function s = jitter2ndx(jit,L,sortflag,Lstart)
%
% using the jitter structure described in jit, creates a jittered
%  spike train using a partition of size L and a randomly chosen starting
%  point for the partition (or Lstart, if it is provided)
%
% jit is a cell array of T different Ntx3 vectors of integers
%   with jit{t}(j,1) <= jit{t}(j,2) < jit{t}(j,3)
% s is a cell array of T different Ntx1 vectors of integers
%
% L is the width of a jitter window (a positive integer)
%
% Lstart is the starting point of the jitter partition
% Lstart defaults to a (discrete uniform) randomly chosen integer in {1,..,L}
% 
% s{t}(j) is a (discrete uniform) random integer in {a,...,b-1}, where
%
%    d = floor((jit{t}(j,2)-Lstart)/L)*L+Lstart;
%    a = max(jit{t}(j,3),d);
%    b = min(jit{t}(j,3),d+L);
%
% note that a <= jit{t}(j,2) < b
%
% sortflag = true (default) sorts s{t}


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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE SPIKE TRAIN INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no cell input
if ~iscell(jit)

    if isempty(jit)
        s = zeros(0,1);
        return
    end
    
    [r,c] = size(jit);
    if c ~= 3
        error('jit must have exactly 3 columns')
    end
    
    s = floor((jit(:,2)-Lstart)./L).*L + Lstart;
    a = max(s,jit(:,1));
    b = min(s+L,jit(:,3));
    s = a + floor(rand(r,1).*(b-a));
    
    if sortflag
        s = sort(s);
    end
    
    return
end % if ~iscell(jit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = cell(size(jit));
T = numel(jit);

for t = 1:T
   
    % extract the data
    jitt = jit{t};
    
    if isempty(jitt)
        s{t} = zeros(0,1);
        continue
    end
    
    [r,c] = size(jitt);
    if c ~= 3
        error(['jit{' num2str(t) '} must have exactly 3 columns'])
    end
    
    st = floor((jitt(:,2)-Lstart)./L).*L + Lstart;
    a = max(st,jitt(:,1));
    b = min(st+L,jitt(:,3));
    st = a + floor(rand(r,1).*(b-a));
    
    if sortflag
        st = sort(st);
    end
    
    s{t} = st;
    
end % loop over trials
