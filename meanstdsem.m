function [meanmat, stdmat, semmat] = meanstdsem(inputmat)
% This function computes the mean, standard deviation, and standard error
% of the mean (SEM) for each column of the input matrix.
% 
% INPUT:
%   inputmat: A matrix where statistics are calculated column-wise.
%
% OUTPUT:
%   meanmat: Row vector containing the mean of each column.
%   stdmat: Row vector containing the standard deviation of each column.
%   semmat: Row vector containing the SEM of each column.

% Calculate the number of columns
numbcol = size(inputmat, 2);     

% Compute the mean of each column, ignoring NaNs
meanmat = nanmean(inputmat, 1);     

% Compute the standard deviation of each column, ignoring NaNs
stdmat = nanstd(inputmat, 0, 1);     

% Count the valid (non-NaN) elements in each column
nmat = sum(~isnan(inputmat), 1);

% Compute the standard error of the mean
semmat = stdmat ./ sqrt(nmat);

end


