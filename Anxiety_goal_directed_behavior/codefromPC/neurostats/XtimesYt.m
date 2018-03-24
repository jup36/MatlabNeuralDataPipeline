function c = XtimesYt(x,y,d,p)
% function c = XtimesYt(x,y,d,p)
%
% x is MxT
% y is NxT
% c = x*(y.')./d and is MxN
%
% d defaults to T (use T-1 for the sample covariance matrix)
%
% If x(:,t) and y(:,t) are the t-th observations of vectors x and y, then
% c(i,j) is the empirical average of x(i)y(j).
%
% If you know that x and y are not binary valued and not mostly zeros, then 
% you might as well use c = (x./d)*(y.').
%
% p is a sparsity (proporition of ones) constraint (default .06) that
% binary inputs must satisfy to use the fast algorithm.


if nargin < 4 || isempty(p)
    p = .06;
elseif numel(p) ~= 1
    error('p must be a scalar')
end

% sizing
[xn,n] = size(x);
[yn,m] = size(y);
if n ~= m, error('x and y muist have the same number of columns'), end

cnumel = xn.*yn;
xnumel = xn.*n;
ynumel = yn.*n;

if nargin < 3 || isempty(d)
    d = 1./n;
else
    d = 1./d;
end

% check for some special cases
if n < 3 || xn*yn < 10 || (x(1)~=(x(1)==1)) || (y(1)~=(y(1)==1)) || sum(x(:)) > p.*xnumel || sum(y(:)) > p.*ynumel
    % probably not worth checking logical
    if xnumel <= ynumel && xnumel <= cnumel
        c = (x.*d)*(y.');
    elseif ynumel <= cnumel
        c = x*((y.*d).');
    else
        c = (double(x)*(y.')).*d;
    end
    return
end

% big and sparse
% quick check for the special logical case
if islogical(x) && islogical(y)

    % initialization
    c = zeros(xn,yn);
    
    % loop through observations
    for k = 1:n
        % extract the data
        xndx = x(:,k);
        yndx = y(:,k);
        % since the data are logical, this syntax gets exactly those
        % elements that should be updated
        c(xndx,yndx) = c(xndx,yndx) + d;
    end

    return

end

% check for logical data
xlog = x == 1;
ylog = y == 1;

if isequal(x,xlog) && isequal(y,ylog)

    % initialization
    c = zeros(xn,yn);
    
    % loop through observations
    for k = 1:n
        % extract the data
        xndx = xlog(:,k);
        yndx = ylog(:,k);
        % since the data are logical, this syntax gets exactly those
        % elements that should be updated
        c(xndx,yndx) = c(xndx,yndx) + d;
    end

    return
end

% do the multiplication
if xnumel <= ynumel && xnumel <= cnumel
    c = (x.*d)*(y.');
elseif ynumel <= cnumel
    c = x*((y.*d).');
else
    c = (double(x)*(y.')).*d;
end
