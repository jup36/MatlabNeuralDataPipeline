function jit = ndx2jitter(s,walls,fixndx,ftype)
%function jit = ndx2jitter(ndx,walls,fixndx,ftype)
%
% creates a data structure that is useful for jittering
%
% ndx is a cell array of spike indices (or a single vector of spike indices)
% jit is a cell array of a Njx3 matrices, where Nj is the number of
% of spikes in ndx{j} (or jit is just a Nx3 matrix when ndx is just a vector)
%
% jit{k}(j,1) is the minimum starting index of any jitter window for spike ndx{k}(j)
% jit{k}(j,2) is the index that every jitter window for spike ndx{k}(j)
%  must contain, namely, ndx{k}(j)
% jit{k}(j,3)-1 is the maximum ending index of any jitter window for spike ndx{k}(j)
%
% walls is a (possibly empty) list of indices that spikes cannot be jittered
%  through (make walls a cell array to use different walls on each trial)
% spike indices that land exactly on a wall are assumed to be to the right of
%  the wall
% a spike is allowed to jitter onto a wall to it's left, but cannot cross it
% a spike is allowed to jitter 1 bin before a wall to it's right, but cannot
%  land on it
% if every ndx is between 1 and M, use walls = [1,M+1] to constrain jittered
%  spikes to stay between 1 and M
%
% fixndx is a cell array of indices into ndx; the default is empty
% the spike indices indicated by fixndx are not jittered
% if fixndx is just a vector of indices, then these same indices are used
%   for all elements of ndx
% +inf can be used to get the last spike index (which is not the largest
%   spike index if s is not sorted)
%
% ftype is a str; the default is 'free'; the options are:
%   'free'      : other spikes can be freely jittered through fixed spikes
%   'wall   '   : fixed spikes are added to the list of walls

%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%

% deal with vector input
nocellF = ~iscell(s);
if nocellF
    s = {s};
end
T = numel(s);

% walls
wallsF = true;
if nargin < 2 || isempty(walls)
    wallsF = false;
end

% convert walls to cell format
if wallsF
    if ~iscell(walls)
        wallstmp = sort(walls);  % saves repeated sorting later
        walls = cell(T,1);
        for t = 1:T
            walls{t} = wallstmp;
        end
    elseif numel(walls) ~= T
        error('s and walls must be similarly sized')
    end
end

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
    elseif numel(fixndx) ~= T
        error('s and fixndx must be similarly sized')
    end
end

% convert inf to last spike time
% convert fixndx to logical
if fixF
    for t = 1:T
        fixndxt = fixndx{t}(:);
        nst = numel(s{t});
        if nst == 0
            fixndx{t} = false(0);
        elseif ~islogical(fixndxt)
            fixndxt(isinf(fixndxt)) = nst;
            f = false(nst,1);
            f(fixndxt) = true;
            fixndx{t} = f;
        elseif nst ~= numel(fixndxt)
            error('sizes of s and fixndx must match when fixndx is logical')
        end
    end
end

% ftype
if nargin < 4 || isempty(ftype)
    ftype = 'free';
end
ftype = lower(ftype);

freeF = strcmp(ftype,'free');
barrierF = strcmp(ftype,'wall');

if ~(freeF || barrierF)
    error('unknown ftype')
end

% treat barrier spikes as walls
if fixF && barrierF && ~wallsF
    wallsF = true;
    walls = cell(T,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JITTER WINDOWS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jit = cell(size(s));

for t = 1:T

    % extract the data from the cells
    st = s{t}(:);
    if ~isequal(st,round(st))
        error('s must be integer valued')
    end
    nt = numel(st);

    if fixF
        fixndxt = fixndx{t}(:);
    end
    if wallsF
        if fixF && barrierF
            wallst = sort([walls{t}(:); st(fixndxt)]);
        else
            wallst = walls{t}(:);
            if ~issorted(wallst)
                wallst = sort(wallst);
            end
        end
        if ~isequal(wallst,round(wallst))
            error('walls must be integer valued')
        end
        nw = numel(wallst);
    end

    % walls
    if wallsF && nw ~= 0

        % sorting
        sortflag = false;
        if ~issorted(st)
            [st,sortndx] = sort(st);
            sortflag = true;
        end

        % initialization
        jitt = [-inf(nt,1) , st , inf(nt,1)];
        
        wallstn = wallst(nw);

        pwval = -inf;
        wj = 1;
        nwval = wallst(1);

        % for each spike, find the walls around it
        for k = 1:nt

            stk = st(k);

            if stk >= nwval % check if the walls need to be incremented

                if stk >= wallstn % check the end
                    
                    jitt(k:nt,1) = nwval;
                    
                    break
                end

                % move walls ahead

                wj = wj + 1;
                nwval = wallst(wj);

                while stk >= nwval

                    wj = wj + 1;
                    nwval = wallst(wj);

                end

                pwval = wallst(wj-1);
            end

            jitt(k,1) = pwval;
            jitt(k,3) = nwval;

        end % loop over spikes

        % unsort
        if sortflag
            st(sortndx) = st;
            jitt(sortndx,:) = jitt;
        end

    else % if wallsF

        % initialization
        jitt = [-inf(nt,1) , st , inf(nt,1)];
    end
        
    % fixed spikes
    if fixF

        % fix the spikes
        jitt(fixndxt,1) = st(fixndxt);
        jitt(fixndxt,3) = st(fixndxt)+1;

    end

    % store 
    jit{t} = jitt;

end % cell loop

% remove artificial cell
if nocellF
    jit = jit{1};
end
