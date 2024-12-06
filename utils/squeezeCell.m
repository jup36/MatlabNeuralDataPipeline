function outputCell = squeezeCell(inputCell)
% squeezeCell removes empty entries from a 1D cell array while preserving
% its original orientation (row or column).
%
% INPUT:
%   inputCell - A 1xN or Nx1 cell array.
%
% OUTPUT:
%   outputCell - A 1x(N-x) or (N-x)x1 cell array with empty entries removed.

    % Validate input
    if ~iscell(inputCell) || (size(inputCell, 1) > 1 && size(inputCell, 2) > 1)
        error('Input must be a 1D cell array (1xN or Nx1).');
    end
    
    % Identify non-empty entries
    nonEmptyMask = ~cellfun(@isempty, inputCell);
    
    % Retain only non-empty entries while preserving orientation
    outputCell = inputCell(nonEmptyMask);
end
