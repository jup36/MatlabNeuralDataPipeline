function [c,md] = jpsth2ccg(J,w,str)
% function [c,N] = jpsth2ccg(J,w,str)
%
% converts a MxM jpsth matrix J to a cross-correlogram c
%
% c is initially a 2M-1 column vector
%
% c(M+k) = mean of elements in J along the kth diagonal, where the
% 0th diagonal is the main diagonal (beginning at (1,1)), the kth diagonal 
% for k < 0 begins at row k+1, i.e., (k+1,1) and the kth diagonal for k > 0
% begins at column k+1, i.e., (1,k+1) 
%
% If w is provided, then only the elements of c corresponding to diagonals 
% -w+1:w-1 are computed and returned, i.e., the central 2w-1 elements of c.
% w defaults to M.
%
% str is an optional string
% the options are
% 
% [] : uses default
% 'full': (default) uses all possible parts of J
% 'rowvalid': the first and last w-1 rows of J are excluded from the
%             computation of c, so that each element of c is a mean over
%             the same number of elements of J
% 'colvalid': the first and last w-1 columns of J are excluded from the
%             computation of c, so that each element of c is a mean over
%             the same number of elements of J
%
% the valid options are only permitted when 2*(w-1) < M
%
% N(k) gives the number of elements of J that were averaged to get c(k)
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(J);

if m ~= n
    error('jpsth2ccg assumes square jpsth')
end

if nargin < 2 || isempty(w)
    w = m;
end

if nargin < 3 || isempty(str)
    str = 'full';
end

str = lower(str);

fullF = strcmp(str,'full');
rowvalidF = strcmp(str,'rowvalid');
colvalidF = strcmp(str,'colvalid');

if ~(fullF || rowvalidF || colvalidF)
    error('unknown str option')
end

if w < 1 || w > m
    error('w must between 1 and M')
end

if ~fullF && 2*(w-1) >= m
    error('w is too large for the valid options')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM (uses DiagSum) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = (1-w:w-1).';

if fullF

    rmin = 1;
    rmax = m;
    cmin = 1;
    cmax = m;
    md = m-abs(d);

elseif rowvalidF

    rmin = w;
    rmax = m-w+1;
    cmin = 1;
    cmax = m;
    md = m-2*w+2;
    
elseif colvalidF

    rmin = 1;
    rmax = m;
    cmin = w;
    cmax = m-w+1;
    md = m-2*w+2;

else
    
    error('unrecognized option')
    
end

c = DiagSum(J,d,rmin,rmax,cmin,cmax) ./ md;

% optional output
if nargout > 1
    
    if numel(md) == 1
        md = repmat(md,size(c));
    end

end