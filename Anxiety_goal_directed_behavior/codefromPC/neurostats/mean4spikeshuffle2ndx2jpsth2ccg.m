function [c,N] = mean4spikeshuffle2ndx2jpsth2ccg(ss,M,L,Lstart,ndy,My,w,cstr,algK,p)
% function [c,N] = mean4spikeshuffle2ndx2jpsth2ccg(ss,M,L,Lstart,ndy,My,w,cstr,algK,p)
%
% computes the expected value of repeated calls to:
%
%   ndx = spikeshuffle2ndx(ss,L,[],Lstart);
%   [c,N] = ndx2jpsth2ccg(ndx,M,ndy,My,'raw',w,cstr,algK);
%
% use [] to get defaults (for example, use Lstart=[] to average over all
%  partitions)
%
% alternative commands that give the same answer are
%
%   p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart);
%   b = ndx2bins(ndy,My);
%   [c,N] = bins2jpsth2ccg(p,b,'raw',w,cstr);
%
% or
%
%   j = mean4spikeshuffle2ndx2jpsth(ss,M,L,Lstart,ndy,My);
%   [c,N] = jpsth2ccg(j,w,cstr);
%
% algK does not affect the answer, but controls the algorithm as in
%  ndx2jpsth2ccg (although here we require (w+L) <= algK*(M+L))
%
% the alternative argument p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING FOR JITTER ARGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = ss{4};

% M check
if numel(M) ~= 1 || M <= 0 || M ~= round(M)
    error('M must be a positive integer scalar')
end

% L check
if numel(L) ~= 1 || L <= 0 || L ~= round(L)
    error('L must be a positive integer scalar')
end

% Lstart check
if nargin < 4 || isempty(Lstart)
    Lstart = [];
elseif numel(Lstart) ~= 1 || Lstart <= 0 || Lstart > L || Lstart ~= round(Lstart)
    error('Lstart must be a scalar integer between 1 and L')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING FOR CCG ARGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(ndy), ndy = {ndy}; end

if numel(ndy) ~= T, error('ndy must have ss{4} elements'), end

if ~isequal(M,My), error('M must equal My because jpsth2ccg assumes square jpsth matrix'), end

if nargin < 7 || isempty(w)
    w = M;
end

if nargin < 8 || isempty(cstr)
    cstr = 'full';
end

cstr = lower(cstr);

fullF = strcmp(cstr,'full');
rowvalidF = strcmp(cstr,'rowvalid');
colvalidF = strcmp(cstr,'colvalid');

if ~(fullF || rowvalidF || colvalidF)
    error('unknown cstr option')
end

if w < 1 || w > M
    error('w must between 1 and M')
end

if ~fullF && 2*(w-1) >= M
    error('w is too large for the valid options')
end

if nargin < 9 || isempty(algK)
    algK = .085;
end
if numel(algK) ~= 1 || algK < 0 || algK > 1
    error('algK must be a scalar in [0,1]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ccg setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w2m1 = 2*w-1;

% decide how to compute the ccg
if w <= algK*M
    compfullF = false;
    % compute only the returned part of the ccg
    c = zeros(w2m1,1);
    % keep all the indices
    keepndx = (1:w2m1).';
else
    compfullF = true;
    % compute the whole raw ccg to avoid if statements inside the loops
    c = zeros(2*M-1,1);
    % only keep the central ndices
    keepndx = (M-w+1:M+w-1).';
end

dg = (1-w:w-1).';

if fullF

    rmin = 1;
    rmax = M;
    cmin = 1;
    cmax = M;
    md = M-abs(dg);
    
elseif rowvalidF

    rmin = w;
    rmax = M-w+1;
    cmin = 1;
    cmax = M;
    md = M-2*w+2;

elseif colvalidF

    rmin = 1;
    rmax = M;
    cmin = w;
    cmax = M-w+1;
    md = M-2*w+2;

else

    error('unrecognized cstr option')

end

if nargout > 1

    if numel(md) == 1
        N = repmat(md,size(dg));
    else
        N = md;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE SHUFFLE MEANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10 || isempty(p)
    p = mean4spikeshuffle2ndx2bins(ss,M,L,Lstart);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE x EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xt = unique(ss{1}(:,1));

% remove excluded data
xt = xt(xt >= rmin & xt <= rmax);
nx = numel(xt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over trials
for t = 1:T

    % extract the data and remove invalid spikes
    yt = ndy{t}(:);
    yt = yt(yt >= cmin & yt <= cmax);
    ny = numel(yt);

    % decide which algorithm to use
    if compfullF

        % compute the full ccg
        % (this avoids if statements inside loops, which is slow)

        yt = yt + M;

        % loop over x spikes
        for j = 1:nx
            
            % get index
            xtj = xt(j);
            % get the contribution
            d = p(xtj,t);
            
            % loop over y spikes
            for k = 1:ny

                % get the entry in the ccg
                ytk = yt(k)-xtj;

                % update the ccg
                c(ytk) = c(ytk) + d;
            end
        end

    else
        
        % compute only the partial ccg

        % assumes sorted spikes
        if ~issorted(yt), yt = sort(yt); end

        yt = yt + w;

        kst = 1;

        % loop over x spikes
        for j = 1:nx
            
            % get index
            xtj = xt(j);
            % get the contribution
            d = p(xtj,t);
        
            % loop over y spikes
            for k = kst:ny
                
                % get the entry in the ccg
                ytk = yt(k)-xtj;
                
                % check if this is a valid index
                if ytk > w2m1
                    break
                elseif ytk < 1
                    kst = kst + 1;
                else
                    % update the ccg
                    c(ytk) = c(ytk) + d;
                end
            end
        end

    end
end

% divide by the number of bins in the jpsth that were summed and only
% keep the relevant indices
c = c(keepndx) ./ (md.*T);
