function [intV,xq] = interpsm(x,v)
xq = x(1):x(end); % new timescale
inthTrj = @(a) interp1(x,a,xq); % interpolation function
vC = mat2cell(v,[ 1 1 1 ], size(v,2)); % convert to cell
intVC = cellfun(@(a) inthTrj(a), vC, 'un', 0); % interpolated hTrj cell
intVCf = cellfun(@(a) sgolayfilt(a,3,201), intVC, 'un', 0); % filtered hTjr cell
intV = cell2mat(intVCf);
end