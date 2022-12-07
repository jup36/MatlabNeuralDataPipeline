function [extV] = extm(x,xq,v,extMethod)
exthTrj = @(a) interp1(x,a,xq,extMethod,'extrap'); % extrapolation function
vC = mat2cell(v,[ 1 1 1 ], size(v,2)); % convert to cell
extVC = cellfun(@(a) exthTrj(a), vC, 'un', 0); % extrapolated hTrj cell
extV = cell2mat(extVC); % extrapolated hTrj
end