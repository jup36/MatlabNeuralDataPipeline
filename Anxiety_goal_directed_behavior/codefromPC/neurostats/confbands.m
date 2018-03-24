function [cbpw,cbsim] = confbands(bsdata,p,meandata,stddata)
% function [cbpw,cbsim] = confbands(bsdata,p,meandata,stddata)
%
% constructs p*100% pointwise (cbpw) and simultaneous (cbsim) confidence
%  bands based on the bootstrapped data (bsdata) and the mean bootstrapped
%  data (meandata)
%
% the kth bootstrap sample is stored in bsdata(:,k)
% cbpw = prctile(bsdata,[(1-p)/2*100 (1-(1-p)/2)*100],2);
%  which is just the pointwise percentiles of bsdata that give p*100% of
%  the data inside the percentiles
%
% cbsim is based on the percentiles of the maximum of the absolute standardized
%  bsdata computed by subtracting meandata from each column and dividing
%  the result by stddata
%
% meandata defaults to all zeros and stddata defaults to all ones

if numel(p) ~= 1 || p < 0 || p > 1
    error('p must be a scalar in [0,1]')
end

[r,c] = size(bsdata);

if nargin < 3 || isempty(meandata)
    meandata = zeros(r,1);
end
if numel(meandata) == 1
    meandata = meandata*ones(r,1);
end

if ~isequal(size(meandata),[r 1])
    error('meandata must the same size as a column of bsdata')
end

if nargin < 4 || isempty(stddata)
    stddata = ones(r,1);
end
if numel(stddata) == 1
    stddata = stddata*ones(r,1);
end
if ~isequal(size(stddata),[r 1])
    error('stddata must the same size as a column of bsdata')
end
if any(stddata <= 0)
    error('stddata must be positive')
end

% compute the pointwise bands
cbpw = prctile(bsdata,[(1-p)/2*100 (1-(1-p)/2)*100],2);

% compute the simultaneous bands
if nargout > 1
    
    % get the p*100th percentile in the empirical distribution of maximum
    % absolute deviations
    m = prctile(max(abs((bsdata - meandata*ones(1,c))./(stddata*ones(1,c)))),p*100);
    
    % create the simultaneous confidence bands
    stdm = m*stddata;
    cbsim = [meandata-stdm, meandata+stdm];
    
end
