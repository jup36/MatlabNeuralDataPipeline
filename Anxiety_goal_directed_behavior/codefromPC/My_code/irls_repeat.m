function [b, iter] = irls_repeat(v,y, xind, b0, offset, lambda)
if nargin<3
    disp('input is not enough');
    return;
else
    if std(v(:,1))~=0
        v = [ones(size(v,1), 1), v];
    end
    vt = v';
    Y = logical(y);
end

if nargin<4 || isempty(b0)
    b0 = rand(size(v,2),1);
end
b_backup = b0; 


if nargin<5 || isempty(offset)
    offset = 0*y;
end
offset_exp = exp(offset);

if nargin<6
    lambda = 0; 
end



% X = v(xind);
logmu_v = v*b0;
logmu = logmu_v(xind)+offset;
mu_v = exp(logmu_v);
mu = exp(logmu);

[numR, numC] = size(v);
temp1 = 0*b0;
for m=1:numC
    temp1(m) = mu'*v(xind,m);
end
temp2 = sum(vt(:, xind(Y)),2);
dmu = temp2-temp1-lambda*b0;
%dmu = xt*(y-mu);
L0 = sum(logmu(Y))-sum(mu)-lambda*(b0')*b0/2;

vij = zeros(numR, numC, numC);
for m=1:numR
    vij(m,:,:) = -v(m,:)'*v(m,:)*mu_v(m);
end
x_ind_count = hist(xind, 1:numR);
for iter=1:100
    h = zeros(numC);
    for m=1:numR
        if x_ind_count(m)>0
            h = h+sum(offset_exp(xind==m))*squeeze(vij(m,:,:));
        end
    end
    h = h-lambda*eye(numC); 
    %     h = xt*bsxfun(@times, mu, x);
    db = h\dmu;
    b = b0-db;
    
    logmu_v = v*b;
    logmu = logmu_v(xind)+offset;
    mu_v = exp(logmu_v);
    mu = exp(logmu);
    
    vij = zeros(numR, numC, numC);
    for m=1:numR
        vij(m,:,:) = -v(m,:)'*v(m,:)*mu_v(m);
    end
    
    for m=1:numC
        temp1(m) = mu'*v(xind,m);
    end
    temp2 = sum(vt(:, xind(Y)),2);
    dmu = temp2-temp1-lambda*b;
    %dmu = xt*(y-mu);
    L = sum(logmu(Y))-sum(mu)-lambda*(b')*b/2;
    
    %     disp([L, L-L0]);
    %     pause;
    if abs(L-L0)<0.001
        return;
    end
    
    b0 = b;
    if sum(isnan(b0))
%        disp('wrong parameters'); 
        b = b_backup; 
        return; 
    end
    L0 = L;
end