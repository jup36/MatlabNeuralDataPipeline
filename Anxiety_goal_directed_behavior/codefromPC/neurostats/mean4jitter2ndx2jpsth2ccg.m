function [c,N] = mean4jitter2ndx2jpsth2ccg(jit,M,L,Lstart,ndy,My,w,cstr,algK)
% function [c,N] = mean4jitter2ndx2jpsth2ccg(jit,M,L,Lstart,ndy,My,w,cstr,algK)
%
% computes the expected value of repeated calls to:
%
%   ndx = jitter2ndx(jit,L,[],Lstart);
%   [c,N] = ndx2jpsth2ccg(ndx,M,ndy,My,'raw',w,cstr,algK);
%
% use [] to get defaults (for example, use Lstart=[] to average over all
%  partitions)
%
% alternative commands that give the same answer are
%
%   p = mean4jitter2ndx2bins(jit,M,L,Lstart);
%   b = ndx2bins(ndy,My);
%   [c,N] = bins2jpsth2ccg(p,b,'raw',w,cstr);
%
% or
%
%   j = mean4jitter2ndx2jpsth(jit,M,L,Lstart,ndy,My);
%   [c,N] = jpsth2ccg(j,w,cstr);
%
% algK does not affect the answer, but controls the algorithm as in
%  ndx2jpsth2ccg (although here we require (w+L) <= algK*(M+L))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING FOR JITTER ARGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no cell input
if ~iscell(jit)
    jit = {jit};
end

T = numel(jit);

% M check
if numel(M) ~= 1 || M <= 0 || M ~= round(M)
    error('M must be a positive integer scalar')
end

% L check
if numel(L) ~= 1 || L <= 0 || L ~= round(L)
    error('L must be a positive integer scalar')
end

% Lstart check
LstartF = true;
if nargin < 4 || isempty(Lstart)
    LstartF = false;
elseif numel(Lstart) ~= 1 || Lstart <= 0 || Lstart > L || Lstart ~= round(Lstart)
    error('Lstart must be a scalar integer between 1 and L')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PROCESSING FOR CCG ARGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(ndy), ndy = {ndy}; end

if numel(ndy) ~= T, error('jit and ndy must have the same number of elements'), end

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

% decide how to compute the ccg
if (w+L) <= algK*(M+L) 
    compfullF = false;
    % compute only the returned part of the ccg
    % (use start and end garbage indices)
    cMaxNdx = 2*w+2*L;
    cMidNdx = w+2*L+1;
    c = zeros(2*w+1+4*L,1);
    % keep all the non-garbage indices
    keepndx = (2*L+2:2*w+2*L).';
else
    compfullF = true;
    % compute the whole raw ccg to avoid if statements inside the loops
    % (use start and end garbage indices)
    cMidNdx = M+1;
    c = zeros(2*M+1,1);
    % only keep the central non-garbage indices
    keepndx = (M-w+2:M+w).';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM LSTART FIXED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the main algorithm works exactly like the 'raw' version of ndx2jpsth2ccg,
% however, the ccg is updated with the mean in a manner analogous to
% mean4jitter2ndx

% there is a fundamental difference in Lstart fixed verses Lstart random

if LstartF

    % loop over trials
    for t = 1:T

        % extract the jitter data
        jitt = jit{t};
        [nt,tmp3] = size(jitt);
        if nt == 0, continue, end
        if tmp3 ~= 3, error(['jit{' num2str(t) '} must be a matrix with 3 columns']), end

        % remove noncontributing spikes
        jitt = jitt(jitt(:,2) >= 2-L & jitt(:,2) <= M+L-1 & jitt(:,1) <= M & jitt(:,3) > 1,:);
        [nt,tmp3] = size(jitt);
        if nt == 0, continue, end        
        
        % get the jitter window for this Lstart
        st = floor((jitt(:,2)-Lstart)./L).*L + Lstart;
        a = max(st,jitt(:,1));
        b = min(st+L,jitt(:,3));

        % get the contribution of each jitter window
        Lt = 1./(T.*(b-a));

        % deal with excluded indices
        a = min(max(a,rmin),rmax+1);
        b = min(max(b,rmin),rmax+1);

        % extract the spike data and remove invalid spikes
        yt = ndy{t}(:);
        yt = yt(yt >= cmin & yt <= cmax);
        ny = numel(yt);
        if ny == 0, continue, end
        % change y's so that a single subtraction gets the ccg ndx
        yt = yt + cMidNdx;

        % decide which algorithm to use
        if compfullF

            % compute the full ccg
            % (this avoids if statements inside loops, which is slow)

            % loop over y spikes
            for k = 1:ny

                ytk = yt(k);

                % loop over x spikes
                for j = 1:nt

                    % get the end entry in the ccg
                    aj = ytk-a(j);
                    % get the start entry in the ccg
                    bj = ytk-b(j);
                    % get the contribution of this spike
                    Ltj = Lt(j);

                    % update the ccg
                    c(bj) = c(bj) - Ltj;
                    c(aj) = c(aj) + Ltj;
                end
            end

        else

            % compute only the partial ccg

            % assumes sorted spikes
            if ~issorted(yt), yt = sort(yt); end
            if ~issorted(b)
                [b,sortndx] = sort(b);
                a = a(sortndx);
                Lt = Lt(sortndx);
            end

            jst = 1;

            % loop over y spikes
            for k = 1:ny

                ytk = yt(k);

                % loop over x spikes
                for j = jst:nt

                    % get the start entry in the ccg
                    bj = ytk-b(j);

                    % check if this is a valid index
                    if bj > cMaxNdx
                        jst = jst + 1;
                    elseif bj < 1
                        break
                    else
                        % get the end entry in the ccg
                        aj = ytk-a(j);
                        % get the contribution of this spike
                        Ltj = Lt(j);

                        % update the ccg
                        c(bj) = c(bj) - Ltj;
                        c(aj) = c(aj) + Ltj;
                    end
                end
            end

        end
    end

    % use cumsum (in reverse) to compute the actual ccg
    c = flipud(cumsum(flipud(c)));

    % divide by the number of bins in the jpsth that were summed and only
    % keep the relevant indices
    c = c(keepndx) ./ md;

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM LSTART RANDOM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Lstart random, we can use cumsum twice to get the triangular
% distribution away from the edges, and then we can loop over all
% possible Lstarts near the edge

Ltri = 1./(T.*L.*L);
Ltri2 = -2.*Ltri;
ctri = zeros(size(c));

% loop over trials
for t = 1:T

    % extract the jitter data
    jitt = jit{t};
    [nt,tmp3] = size(jitt);
    if nt == 0, continue, end
    if tmp3 ~= 3, error(['jit{' num2str(t) '} must be a matrix with 3 columns']), end

    % remove noncontributing spikes
    jitt = jitt(jitt(:,2) >= 2-L & jitt(:,2) <= M+L-1 & jitt(:,1) <= M & jitt(:,3) > 1,:);
    [nt,tmp3] = size(jitt);
    if nt == 0, continue, end

    % extract the spike data and remove invalid spikes
    yt = ndy{t}(:);
    yt = yt(yt >= cmin & yt <= cmax);
    ny = numel(yt);
    if ny == 0, continue, end
    
    % change y's so that a single subtraction gets the ccg ndx
    yt = yt + cMidNdx;

    % sort if necessary
    if ~compfullF
        % assumes sorted spikes
        if ~issorted(yt), yt = sort(yt); end
        if ~issorted(jitt(:,2))
            jitt = sortrows(jitt,2);
        end
    end
    
    % find those jitter windows that give triangular distirbutions    
    j1 = jitt(:,2)-(L-1);
    j3 = j1+(2*L-1);
    triind = (jitt(:,1) <= j1 & j3 <= jitt(:,3) & j1 >= rmin & j3 <= rmax+1); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the triangular windows
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a = j1(triind);
    ab = a + L;
    b = ab + L;
    nt = numel(a);
    
    if nt > 0
    
        % decide which algorithm to use
        if compfullF

            % compute the full ccg
            % (this avoids if statements inside loops, which is slow)

            % loop over y spikes
            for k = 1:ny

                ytk = yt(k);

                % loop over x spikes
                for j = 1:nt

                    % get the end entry in the ccg
                    aj = ytk-a(j);
                    % get the middle entry in the ccg
                    abj = ytk-ab(j);
                    % get the start entry in the ccg
                    bj = ytk-b(j);

                    % update the ccg
                    ctri(bj) = ctri(bj) + Ltri;
                    ctri(abj) = ctri(abj) + Ltri2;
                    ctri(aj) = ctri(aj) + Ltri;
                end
            end

        else

            % compute only the partial ccg

            % assumes sorted spikes
            if ~issorted(b)
                [b,sortndx] = sort(b);
                a = a(sortndx);
                ab = ab(sortndx);
            end

            jst = 1;

            % loop over y spikes
            for k = 1:ny

                ytk = yt(k);

                % loop over x spikes
                for j = jst:nt

                    % get the start entry in the ccg
                    bj = ytk-b(j);

                    % check if this is a valid index
                    if bj > cMaxNdx
                        jst = jst + 1;
                    elseif bj < 1
                        break
                    else
                        % get the end entry in the ccg
                        aj = ytk-a(j);
                        % get the middle entry in the ccg
                        abj = ytk-ab(j);

                        % update the ccg
                        ctri(bj) = ctri(bj) + Ltri;
                        ctri(abj) = ctri(abj) + Ltri2;
                        ctri(aj) = ctri(aj) + Ltri;
                    end
                end
            end

        end % compfullF
    end % if nt > 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the edge windows
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % copy the data
    
    triind = ~triind;
    j1 = jitt(triind,1);
    j2 = jitt(triind,2);
    j3 = jitt(triind,3);
    
    nt = numel(j1);
    
    if nt > 0
    
        for Ls = 1:L

            % get the jitter window for this Lstart
            st = floor((j2-Ls)./L).*L + Ls;
            a = max(st,j1);
            b = min(st+L,j3);

            % get the contribution of each jitter window
            Lt = 1./((L.*T).*(b-a));

            % deal with excluded indices
            a = min(max(a,rmin),rmax+1);
            b = min(max(b,rmin),rmax+1);


            % decide which algorithm to use
            if compfullF

                % compute the full ccg
                % (this avoids if statements inside loops, which is slow)

                % loop over y spikes
                for k = 1:ny

                    ytk = yt(k);

                    % loop over x spikes
                    for j = 1:nt

                        % get the end entry in the ccg
                        aj = ytk-a(j);
                        % get the start entry in the ccg
                        bj = ytk-b(j);
                        % get the contribution of this spike
                        Ltj = Lt(j);

                        % update the ccg
                        c(bj) = c(bj) - Ltj;
                        c(aj) = c(aj) + Ltj;
                    end
                end

            else

                % compute only the partial ccg

                % assumes sorted spikes
                if ~issorted(yt), yt = sort(yt); end
                if ~issorted(b)
                    [b,sortndx] = sort(b);
                    a = a(sortndx);
                    Lt = Lt(sortndx);
                end

                jst = 1;

                % loop over y spikes
                for k = 1:ny

                    ytk = yt(k);

                    % loop over x spikes
                    for j = jst:nt

                        % get the start entry in the ccg
                        bj = ytk-b(j);

                        % check if this is a valid index
                        if bj > cMaxNdx
                            jst = jst + 1;
                        elseif bj < 1
                            break
                        else
                            % get the end entry in the ccg
                            aj = ytk-a(j);
                            % get the contribution of this spike
                            Ltj = Lt(j);

                            % update the ccg
                            c(bj) = c(bj) - Ltj;
                            c(aj) = c(aj) + Ltj;
                        end
                    end
                end

            end % compfullF
        end % Ls loop
    end % if nt > 0
end % trial loop

% use cumsum (in reverse) to compute the actual ccg
% use cumsum twice for the triangular part
% add the results together
c = flipud(cumsum(flipud(c)+cumsum(flipud(ctri))));

% divide by the number of bins in the jpsth that were summed and only
% keep the relevant indices
c = c(keepndx) ./ md;
