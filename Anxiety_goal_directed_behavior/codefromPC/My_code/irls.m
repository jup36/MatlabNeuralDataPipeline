function [b, iter] = irls(x,y, b0, offset, lambda)
%% fit GLM model with IRLS method
if nargin<2
    disp('input is not enough');
    return;
else
    if std(x(:,1))~=0
        x = [ones(size(x,1), 1), x];
    end
    xt = x';
    Y = logical(y);
end
num_p = size(x,2); 

if nargin<3 || isempty(b0)
    b0 = rand(size(x,2),1);
end

if nargin<4 || isempty(offset)
    offset = 0*y;
end

if nargin<5 || isempty(lambda)
    lambda = 0; 
end

logmu = x*b0 +offset;
mu = exp(logmu);
dmu = xt*(y-mu)-lambda*b0;
L0 = sum(logmu(Y))-sum(mu)-lambda*(b0'*b0)/2;

for iter=1:100
    h = -xt*bsxfun(@times, x, mu)-lambda*eye(num_p);
    db = h\dmu;
    b = b0-db;
    
    logmu = x*b +offset;
    mu = exp(logmu);
    dmu = xt*(y-mu)-lambda*b;
    L = sum(logmu(Y))-sum(mu)-lambda*(b'*b)/2;
    
    %     disp([L, L-L0]);
    %     pause;
    if abs(L-L0)<0.00001
        return;
    end
    
    b0 = b;
    L0 = L;
end