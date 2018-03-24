function s = ndx2psthSD(ndx,Mx,str,p)
% function s = ndx2psthSD(ndx,Mx,str,xpsth)
%
% s = bins2psthSD(ndx2bins(ndx,Mx),str,...);

T = numel(ndx);

s = zeros(Mx,1);

% deal with a single spike train
if ~iscell(ndx) || T == 1    
    return
end

% input processing
if nargin < 3 || isempty(str)
    str = 'raw'; 
end

if nargin < 4 || isempty(p)
    p = ndx2psth(ndx,Mx); 
end

str = lower(str);

rawF = strcmp(str,'raw');
ubF = strcmp(str,'ub');
rawNSF = strcmp(str,'rawns');
ubNSF = strcmp(str,'ubns');

% normal cell routine
d = sqrt(1./T);

for t = 1:T
    
    st = zeros(Mx,1);
    
    ndxt = ndx{t}(:);
    ndxt = ndxt(ndxt >= 1 & ndxt <= Mx);
    numndxt = numel(ndxt);
    
    for k = 1:numndxt
        
        ndxtk = ndxt(k);
        
        st(ndxtk) = st(ndxtk) + d;
        
    end
    
    % this line takes advantage of the fact the matlab ignores multiple
    % indices
    s(ndxt) = s(ndxt) + st(ndxt).^2;
    
end

% subtract the mean and take the square root
if rawF
    
    s = sqrt(s - p.^2);
    return
end

if ubF
    
    s = sqrt((T./(T-1)).*(s - p.^2));
    return
end

error('unrecognized str option')
