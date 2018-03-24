function [ycdf, xquantile] = ksplot(tsp, p, t0t1)
%compute ks plot
if iscell(tsp)
    N = length(tsp);
    if size(p,2)==1
        p = reshape(p, [], N);
    end
else
    N = 1;
    p = reshape(p, [], 1);
end



if exist('t0t1', 'var')
    t0 = t0t1(1);
    t1 = t0t1(2);
else
    t0 = 1/1000;
    t1 = length(p)/1000;
end

xquantile = [];
for m=1:N
    ts = tsp{m};
    ts(ts<t0) = [];
    ts(ts>t1) = [];
    z = zk(ts, p(:,m));
    xquantile = [xquantile; z]; %#ok<AGROW>
end

xquantile = sort(xquantile);
T = length(xquantile);
if T>1
    ycdf = (0.5:(T-0.5))/T;
else
    ycdf = [];
end

function z = zk(ts, p)
%transform the rescaled t
Fs = 1000;
T = length(p);
ts = unique(sort(ceil(ts*Fs)));
ts(isnan(ts)) = [];
ts(ts<0) = [];
ts(ts>T) = [];
if isnan(ts)
    z = [];
    return;
end

temp = cumsum(p);
rescale_t = temp(ts);
temp = reshape(rescale_t, [], 1);
if isempty(temp)
    z =[];
    return;
else
    %tauk = diff(temp);
    tauk = [temp(1); diff(temp)];
end

z = 1-exp(-tauk);