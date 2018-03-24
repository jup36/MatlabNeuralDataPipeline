function yhat = xbeta(x, beta, xind)
%X*beta 
if size(x,2)==length(beta)
    y = x*beta; 
elseif size(x,2)==length(beta)-1
    y = x*beta(2:end)+beta(1); 
else
    disp('size of x does not match with beta!'); 
    yhat = []; 
    return; 
end

if exist('xind', 'var')
    yhat = y(xind); 
else
    yhat = y; 
end

