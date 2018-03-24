function s = bins2psthSD(x,str,p)
% function s = bins2psthSD(x,str,xpsth)
%
% computes the empirical standard deviation of the peri-stimulus time histogram of 
% the binned spike train x, where the kth trial is in the kth column of x
%
% s = std(x,1,2) is a column vector
%
% s(k) = empirical standard deviation of the firing rate of x at time k
%
% xpsth = bins2psth(x) is optional (and is used to simplify calculations)
% (Note that if xpsth is provided, but xpsth ~= bins2psth(x), then the
% results will be incorrect.)
%
% str is string that specifies what type of standard deviation to compute
% the options are
%
% [] : uses default
% 'raw' : (default) normalizes by n (number of columns of x)
% 'ub' : unbiased sample standard deviation (normalizes by n-1)
% 'rawNS' : nonstationary correction for raw
% 'ubNS' : nonstationary correction for ub

if nargin < 2 || isempty(str)
    str = 'raw'; 
end

if nargin < 3 || isempty(p)
    p = bins2psth(x); 
end

str = lower(str);

rawF = strcmp(str,'raw');
ubF = strcmp(str,'ub');
rawNSF = strcmp(str,'rawns');
ubNSF = strcmp(str,'ubns');

n = size(x,2);

% simpler formula for logical x
if islogical(x)
    
    if rawF
        
        s = sqrt(p.*(1-p));
        return
    end
        
    if ubF
        
        s = sqrt((n./(n-1)).*(p.*(1-p)));
        return
        
    end
        
end

% standard (shortcut) formula
if rawF

    s = sqrt(sum(x.^2,2)./n - p.^2);
    return
    
end

if ubF
    
    s = sqrt((n./(n-1)).*(sum(x.^2,2)./n - p.^2));
    return
    
end
