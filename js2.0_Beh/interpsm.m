function [intV,xq] = interpsm(x,v)
xq = x(1):x(end); % new timescale
inthTrj = @(a) interp1(x,a,xq); % interpolation function
vC = mat2cell(v,[ 1 1 1 ], size(v,2)); % convert to cell
intVC = cellfun(@(a) inthTrj(a), vC, 'un', 0); % interpolated hTrj cell

lengthIntVC = unique(cellfun(@length, intVC)); 

% to use sg filter get sgfiltFramelen
if lengthIntVC >= 201
    sgfiltFramelen = 101;
elseif lengthIntVC < 201
    if mod(lengthIntVC,2)==0
        sgfiltFramelen = lengthIntVC-1; % the frame length for sg filter needs to be an odd number
    else
        sgfiltFramelen = lengthIntVC;
    end
end
intVCf = cellfun(@(a) sgolayfilt(a,3,sgfiltFramelen), intVC, 'un', 0); % filtered hTjr cell
intV = cell2mat(intVCf);
end