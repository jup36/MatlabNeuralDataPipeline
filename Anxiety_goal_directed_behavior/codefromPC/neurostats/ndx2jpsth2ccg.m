function [c,N] = ndx2jpsth2ccg(ndx,M,ndy,My,str,w,cstr,algK,xpsth,ypsth,xpsthSD,ypsthSD)
% function [c,N] = ndx2jpsth2ccg(ndx,Mx,ndy,My,str,w,cstr,algK,xpsth,ypsth,xpsthSD,ypsthSD)
%
% J = ndx2jpsth(ndx,Mx,ndy,My,str,xpsth,ypsth,xpsthSD,ypsthSD);
% [c,N] = jpsth2ccg(J,w,cstr);
%
% the arguments following cstr are optional
%
% This can provide a shortcut for computing ccg's by avoiding explicit 
% construction of the jpsth.
%
% algK (default = .085) is a constant between 0 and 1 such that whenever
% w <= algK*M a different algorithm is used.  This can improve performance
% when w is much smaller than M.
% 

if ~iscell(ndx), ndx = {ndx}; end
if ~iscell(ndy), ndy = {ndy}; end

T = numel(ndx);

if numel(ndy) ~= T, error('ndx and ndy must have the same number of elements'), end

if ~isequal(M,My), error('Mx must equal My because jpsth2ccg assumes square jpsth matrix'), end

%%%%%%%%%%%%%%%%%%%%%%
% jpsth preprocessing
%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5 || isempty(str)
    str = 'raw';
end

str = lower(str);

rawF = strcmp(str,'raw');
covF = strcmp(str,'cov');
covubF = strcmp(str,'covub');
corF = strcmp(str,'cor');
corubF = strcmp(str,'corub');
covNSF = strcmp(str,'covns');
covubNSF = strcmp(str,'covubns');
corNSF = strcmp(str,'corns');
corubNSF = strcmp(str,'corubns');

if covNSF || covubNSF || corNSF || corubNSF
    error('option unimplemented')
end
if ~(rawF || covF || covubF || corF || corubF || covNSF || covubNSF || corNSF || corubNSF)
    error('unrecognized option')
end

%%%%%%%%%%%%%%%%%%%%%%
% ccg preprocessing
%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6 || isempty(w)
    w = M;
end

if nargin < 7 || isempty(cstr)
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

if nargin < 8 || isempty(algK)
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
% raw init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rawF
    
    d = 1./T;
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % psths
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 9 || isempty(xpsth)

        xpsth = ndx2psth(ndx,M);

    end

    if nargin < 10 || isempty(ypsth)

        ypsth = ndx2psth(ndy,My);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cov init
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if covF
        
        d = 1./T;
        
        c(keepndx) = DiagSumX(-xpsth,ypsth,dg,rmin,rmax,cmin,cmax);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % covub init
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if covubF
        
        d = 1./(T-1);
        
        c(keepndx) = DiagSumX(xpsth.*(-T./(T-1)),ypsth,dg,rmin,rmax,cmin,cmax);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw, cov, covub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rawF || covF || covubF

    % loop over trials
    for t = 1:T

        % extract the data and remove invalid spikes
        xt = ndx{t}(:);
        xt = xt(xt >= rmin & xt <= rmax);
        nx = numel(xt);

        yt = ndy{t}(:);
        yt = yt(yt >= cmin & yt <= cmax);
        ny = numel(yt);

        % decide which algorithm to use
        if compfullF
            
            % compute the full ccg
            % (this avoids if statements inside loops, which is slow)
        
            yt = yt + M;

            % loop over y spikes
            for k = 1:ny

                ytk = yt(k);

                % loop over x spikes
                for j = 1:nx

                    % get the entry in the ccg
                    xtj = ytk-xt(j);

                    % update the ccg
                    c(xtj) = c(xtj) + d;
                end
            end
            
        else
            
            % compute only the partial ccg
            
            % assumes sorted spikes
            if ~issorted(xt), xt = sort(xt); end
            if ~issorted(yt), yt = sort(yt); end
            
            yt = yt + w;
            
            jst = 1;
            
            % loop over y spikes
            for k = 1:ny

                ytk = yt(k);

                % loop over x spikes
                for j = jst:nx

                    % get the entry in the ccg
                    xtj = ytk-xt(j);
                    
                    % check if this is a valid index
                    if xtj > w2m1
                        jst = jst + 1;
                    elseif xtj < 1
                        break
                    else
                        % update the ccg
                        c(xtj) = c(xtj) + d;
                    end
                    
                end
            end            
            
        end
    end

    % divide by the number of bins in the jpsth that were summed and only
    % keep the relevant indices
    c = c(keepndx) ./ md;

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psth stdard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 11 || isempty(xpsthSD)
    
    xpsthSD = [];
    
    if corF
        xpsthSD = ndx2psthSD(ndx,M,'raw',xpsth);
    end
    
    if corubF
        xpsthSD = ndx2psthSD(ndx,M,'ub',xpsth);
    end
        
end

if nargin < 12 || isempty(ypsthSD)

    ypsthSD = [];
    
    if corF
        ypsthSD = ndx2psthSD(ndy,My,'raw',ypsth);
    end
    
    if corubF
        ypsthSD = ndx2psthSD(ndy,My,'ub',ypsth);
    end

end

% deal with divide by zero
xpsthSD = xpsthSD + sqrt(eps(xpsthSD));
ypsthSD = ypsthSD + sqrt(eps(ypsthSD));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cor init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corF

    dx = 1./(T.*xpsthSD);
    dy = 1./ypsthSD;

    c(keepndx) = DiagSumX(-xpsth./xpsthSD,(ypsth./ypsthSD).',dg,rmin,rmax,cmin,cmax);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corub init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corubF

    dx = 1./((T-1).*xpsthSD);
    dy = 1./ypsthSD;

    c(keepndx) = DiagSumX((xpsth.*(-T./(T-1)))./xpsthSD,(ypsth./ypsthSD).',dg,rmin,rmax,cmin,cmax);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cor, corub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if corF || corubF
    
    % loop over trials
    for t = 1:T

        % extract the data and remove invalid spikes
        xt = ndx{t}(:);
        xt = xt(xt >= rmin & xt <= rmax);
        nx = numel(xt);

        yt = ndy{t}(:);
        yt = yt(yt >= cmin & yt <= cmax);
        ny = numel(yt);
        
        % decide which algorithm to use
        if compfullF
            
            % compute the full ccg
            % (this avoids if statements inside loops, which is slow)

            % loop over y spikes
            for k = 1:ny

                % get the y entry and its normalization
                ytk = yt(k);
                dyk = dy(ytk);
                % get it ready for the ccg entry
                ytk = ytk + M;

                % loop over x spikes
                for j = 1:nx

                    % get the x entry and its normalization
                    xtj = xt(j);
                    d = dx(xtj).*dyk;

                    % get the entry in the ccg
                    xtj = ytk-xtj;

                    % update the ccg
                    c(xtj) = c(xtj) + d;
                end
            end
            
        else
            
            % compute only the partial ccg
            
            % assumes sorted spikes
            if ~issorted(xt), xt = sort(xt); end
            if ~issorted(yt), yt = sort(yt); end
            
            jst = 1;
            
            % loop over y spikes
            for k = 1:ny

                % get the y entry and its normalization
                ytk = yt(k);
                dyk = dy(ytk);
                % get it ready for the ccg entry
                ytk = ytk + w;

                % loop over x spikes
                for j = jst:nx

                    % get the x entry and its index into the ccg
                    xtjtrue = xt(j);
                    xtj = ytk-xtjtrue;
                    
                     % check if this is a valid index
                    if xtj > w2m1
                        jst = jst + 1;
                    elseif xtj < 1
                        break
                    else
                        % get the normalization
                        d = dx(xtjtrue).*dyk;
                        % update the ccg
                        c(xtj) = c(xtj) + d;
                    end
                end
            end            
            
        end
    end

    % divide by the number of bins in the jpsth that were summed and only
    % keep the relevant indices
    c = c(keepndx) ./ md;

    return

end

error('algorithm error: unreachable code')
