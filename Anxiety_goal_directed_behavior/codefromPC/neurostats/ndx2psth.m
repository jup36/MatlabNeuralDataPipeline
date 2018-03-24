function p = ndx2psth(ndx,M)
% function p = ndx2psth(ndx,M)
%
% p = bins2psth(ndx2bins(ndx,M));

p = zeros(M,1);

% deal with a single spike train to avoid cast to cell
if ~iscell(ndx)
    
    ndx = ndx(ndx >= 1 & ndx <= M);
    numndx = numel(ndx);
    
    for t = 1:numndx
        ndxt = ndx(t);
        p(ndxt) = p(ndxt) + 1;
    end
    
    return

end

% normal cell routine
T = numel(ndx);

d = 1./ T;

for t = 1:T
    
    ndxt = ndx{t}(:);
    ndxt = ndxt(ndxt >= 1 & ndxt <= M);
    numndxt = numel(ndxt);
    
    for k = 1:numndxt
        
        ndxtk = ndxt(k);
        
        p(ndxtk) = p(ndxtk) + d;
        
    end
    
end
