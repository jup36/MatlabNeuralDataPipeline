function [b, mu, iter] = irls_long(x,y, b0, offset, lambda)
if nargin<2
    disp('input is not enough');
    return;
else
    num_trial = x.ntrial;
    xbasis = x.fun;
    num_basis = length(xbasis);
    xind = x.ind;
    Y = logical(y);
    T = length(Y)/num_trial;
end

if nargin<3 || isempty(b0)
    b0 = rand(num_basis+1,1);
end

if nargin<4 || isempty(offset)
    offset = 0*y;
end

if nargin<5
    lambda = 0; 
end

logmu = repmat(largeX(xbasis,xind, b0, T), num_trial, 1)+offset;
mu = exp(logmu);
summu = sum(reshape(mu, [], num_trial),2);
dmu = largeXdmu(xbasis,xind, y-mu, num_trial)-lambda*b0;
L0 = sum(logmu(Y))-sum(mu)-lambda*(b0')*(b0)/2;

for iter=1:100
    h = largeXh(xbasis, xind, summu)-lambda*eye(num_basis+1);
    %     h = xt*bsxfun(@times, mu, x);
    db = h\dmu;
    b = b0+db;
    
    logmu = repmat(largeX(xbasis,xind, b, T), num_trial, 1)+offset;
    mu = exp(logmu);
    summu = sum(reshape(mu, [], num_trial),2);
    dmu = largeXdmu(xbasis,xind, y-mu, num_trial)-lambda*b;
    L = sum(logmu(Y))-sum(mu)-lambda*(b')*b/2;
    
    %     disp([L, L-L0]);
    %     pause;
    if (L-L0)<10^(-5)
        mu = mu(1:T); 
        return;
    end
    
    b0 = b;
    L0 = L;
end

function xb = largeX(x,ind,b, T)
if nargin<3
    disp('input is not enough');
    return;
end

if ~exist('T', 'var')
    T = max(ind(:,2));
end

xb = ones(T,1)*b(1);
for m=1:length(x)
    ind0 = ind(m,1);
    ind1 = ind(m,2);
    xb(ind0:ind1) = xb(ind0:ind1)+x{m}*b(m+1);
end

function dmu = largeXdmu(x, ind, y, num_trial)
if nargin<4
    disp('input is not enough');
    return;
end

num_basis = length(x);
dmu = zeros(num_basis+1,1);
y = reshape(sum(reshape(y, [], num_trial),2), 1, []);

dmu(1) = sum(y);
for m=1:num_basis
    ind0 = ind(m,1);
    ind1 = ind(m,2);
    dmu(m+1) = y(ind0:ind1)*x{m};
end


function h = largeXh(xbasis, xind, mu)
if nargin<3
    disp('input is not enough');
    return;
end
num_basis = length(xbasis);
h = zeros(num_basis+1);
h(1,1) = sum(mu); 
for m=2:(num_basis+1)
    x1 = xbasis{m-1};
    ind1_0 = xind(m-1,1);
    ind1_1 = xind(m-1,2);
    h(1,m) = sum(mu(ind1_0:ind1_1).*x1);
    h(m,1) = h(1,m);
    for n= m:num_basis+1
        x2 = xbasis{n-1};
        ind2_0 = xind(n-1,1);
        %         ind2_1 = xind(n-1,1);
        
        if ind2_0>ind1_1
            break;
        end
        tmp_mu = mu(ind2_0:ind1_1);
        dind = ind1_1-ind2_0;
        x1 = x1((end-dind):end);
        x2 = x2(1:(dind+1));
        h(m,n) = sum(x1.*x2.*tmp_mu);
        h(n,m) = h(m,n);
    end
end

















