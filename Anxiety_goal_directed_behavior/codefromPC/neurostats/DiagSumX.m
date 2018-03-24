function a = DiagSumX(x,y,d,rmin,rmax,cmin,cmax)
% function a = DiagSumX(x,y,d,rmin,rmax,cmin,cmax)
%
% x is a R x 1 vector
% y is a C x 1 vector
%
% a = DiagSum(x*y.',d,rmin,rmax,cmin,cmax);
%
% this usage can be faster than the above command because the explicit
% construction of x*y.' is avoided

%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%

% at the end, d should be a reverse sorted column vector, with d1 = max(d)
% and dn = min(d) and n = numel(d)
%
% a should be a matrix the same size as d
%
% if d is a list of consecutive integers, then cFlag should be thrown

R = numel(x);
C = numel(y);

sortflag = false;
cFlag = false;

if nargin < 7 || isempty(cmax)
    cmax = C;
end
if nargin < 6 || isempty(cmin)
    cmin = 1;
end
if nargin < 5 || isempty(rmax)
    rmax = R;
end
if nargin < 4 || isempty(rmin)
    rmin = 1;
end

if rmin > rmax || cmin > cmax || rmin < 1 || rmax > R || cmin < 1 || cmax > C
    error('optional sizing arguments are invalid')
end

if nargin < 3 || isempty(d)

    d1 = C-1;
    dn = 1-R;
    n = C+R-1;
    cFlag = true;

    a = zeros(n,1);

else

    a = zeros(size(d));

    d = d(:);
    n = numel(d);
    dn = d(1);
    d1 = d(n);

    if isequal(d,(dn:d1).')
        cFlag = true;
    else

        if ~issorted(d)
            [d,sortndx] = sort(d);
            sortflag = true;
        end

        d1 = d(n);
        dn = d(1);

        if isequal(d,(dn:d1).')
            cFlag = true;
        else
            % reverse order of d so that it is decreasing
            d(n:-1:1) = d;
        end
    end

    if d1 >= C || dn <= -R
        error('invalid diagonal index d')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%

% min and max columns that contribute to the answer
kstart = max(dn+1,cmin);
kend = min(R+d1,cmax);

% special algorithm for consecutive d's
if cFlag

    d11 = d1+1;
    d1min = d11+rmin;
    d1max = d11+rmax;

    % loop over columns in M
    for k = kstart:kend

        jstart = max(1,d1min-k);
        jend = min(n,d1max-k);

        kd1 = k - d11;
        yk = y(k);

        % loop over items in d
        for j = jstart:jend

            kdj = kd1+j;

            a(j) = a(j) + x(kdj).*yk;
        end

    end

else % general algorithm

    % jstart and jend tell which element of d to consider at any given k
    jstart = n+1;
    jend = n;
    rmin1 = rmin-1;

    % do the remaining columns in M
    for k = kstart:kend

        % figure out jstart and jend
        % move backwards from start to see if more should be included
        for j = jstart-1:-1:1
            % k - d(j) should be decreasing
            if k - d(j) > rmin1
                jstart = j;
            else
                break
            end
        end
        % move backwards from jend to see if more should be included
        for j = jend:-1:1
            % k - d(j) should be decreasing
            if k - d(j) > rmax
                jend = j-1;
            else
                break
            end
        end

        yk = y(k);

        % loop over items in d
        for j = jstart:jend

            kdj = k - d(j);

            a(j) = a(j) + x(kdj).*yk;
        end
    end

end

% flip back
a(n:-1:1) = a;

% unsort
if sortflag
    a(sortndx) = a;
end

