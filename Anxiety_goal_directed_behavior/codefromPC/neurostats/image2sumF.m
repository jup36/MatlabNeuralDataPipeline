function M = image2sumF(J,fsize,shape)
% function M = image2sumF(J,fsize,shape)
%
% J is r x c matrix (presumably from J = image2sum(I);)
%
% fsize = [fr fc] is the size of a filter
% M is a (r+fr-1) x (c+fc-1) matrix with the property that
%
% M(m,n) = J(m',n')-J(m-fr,n')-J(m',n-fc)+J(m-fr,n-fc)
% where m' = min(m,r) and n' = min(n,r)
%
% if m-fr < 1 or n-fc < 1, then the respective terms are zero
%
% See image2sum.m for an explanation of why this is useful.
%
% currently does not work for any(fsize < size(J))
%
% The optional argument shape (default = 'valid') works just like
% filter2.

if nargin < 3
    shape = 0;
elseif ischar(shape)
    switch lower(shape)
        case 'valid'
            shape = 0;
        case 'same'
            shape = 1;
        case 'full'
            shape = 2;
        otherwise
            error('bad shape')
    end
end

[r,c] = size(J);
fr = fsize(1);
fc = fsize(2);

Mr = r+fr-1;
Mc = c+fc-1;

M = zeros(Mr,Mc);

% M(x,y) = J(x',y')-J(x-fr,y')-J(x',y-fc)+J(x-fr,y-fc)

% only compute what is necessary
if shape == 0
    xs = fr;
    xe = r;
    ys = fc;
    ye = c;
elseif shape == 1
    xs = ceil((fr-1)/2)+1;
    xe = xs+r-1;
    ys = ceil((fc-1)/2)+1;
    ye = ys+c-1;
elseif shape == 2
    xs = 1;
    xe = Mr;
    ys = 1;
    ye = Mc;
else
    error('bad shape')
end

M(xs:fr,ys:fc) = J(xs:fr,ys:fc);

for y = ys:fc
    % y-fc terms are gone
    for x = fr+1:r
        M(x,y) = J(x,y)-J(x-fr,y);
    end
    Jry = J(r,y);
    for x = r+1:xe
        % use x' = r
        M(x,y) = Jry-J(x-fr,y);
    end
end

for y = fc+1:c
    for x = xs:fr
        % x - fr terms are gone
        M(x,y) = J(x,y)-J(x,y-fc);
    end
    for x = fr+1:r
        M(x,y) = J(x,y)-J(x-fr,y)-J(x,y-fc)+J(x-fr,y-fc);
    end
    Jry = J(r,y)-J(r,y-fc);
    for x = r+1:xe
        % use x' = r
        M(x,y) = Jry-J(x-fr,y)+J(x-fr,y-fc);
    end
end

for y = c+1:ye
    % use y' = c
    for x = xs:fr
        % x - fr terms are gone
        M(x,y) = J(x,c)-J(x,y-fc);
    end
    for x = fr+1:r
        M(x,y) = J(x,c)-J(x-fr,c)-J(x,y-fc)+J(x-fr,y-fc);
    end
    Jry = J(r,c)-J(r,y-fc);
    for x = r+1:xe
        % use x' = r
        M(x,y) = Jry-J(x-fr,c)+J(x-fr,y-fc);
    end
end

% resize
if shape ~= 2
    M = M(xs:xe,ys:ye);
end

