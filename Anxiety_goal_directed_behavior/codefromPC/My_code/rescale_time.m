function tauk = rescale_time(tsp, p, t0t1)
%% rescale time

if iscell(tsp)
    N = length(tsp);
    p = reshape(p,[], N);
else
    N = 1;
    temp = tsp; 
    tsp = cell(1); 
    tsp{1} = temp; 
    p = reshape(p, [], 1);
end

t0 = t0t1(1); 
t1 = t0t1(2); 
tauk = cell(N, 1); 
for m=1:N
    ts = tsp{m}; 
    ts(ts<=t0) = []; 
    ts(ts>t1) = []; 
    ind = ceil((ts-t0)*1000); 
    temp = cumsum(p(:, m)); 
    L = temp(ind); 
    if isempty(L)
        tauk{m} = []; 
    else
        tauk{m} = [L(1); diff(L)];
    end
end

%tauk = cell2mat(tauk); 
