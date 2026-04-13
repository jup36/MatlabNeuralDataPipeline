function mask = match_time_indices(tFullX, tPartX, tol)
%MATCH_TIME_INDICES Return logical mask for matching timestamps
%
%   mask = MATCH_TIME_INDICES(tFullX, tPartX, tol)
%
% INPUTS:
%   tFullX : reference time vector (1 x N)
%   tPartX : subset of time points to match (1 x M)
%   tol    : optional tolerance for matching (default: 1e-10)
%
% OUTPUT:
%   mask   : logical array (1 x N), true where tFullX matches any tPartX
%
% Example:
%   tFullX = -0.9:0.01:5;
%   tPartX = [0 0.5 1];
%   mask = match_time_indices(tFullX, tPartX);
%

    if nargin < 3
        tol = 1e-10; % safeguard against floating-point precision issues
    end
    
    mask = false(size(tFullX));
    for i = 1:numel(tPartX)
        mask = mask | abs(tFullX - tPartX(i)) < tol;
    end
end
