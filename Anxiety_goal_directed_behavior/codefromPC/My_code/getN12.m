function [N12, tsp12] = getN12(tsp1, tsp2, t0t1, deltat)
%% compute the number of synchronized spikes
t0 = t0t1(1);
t1 = t0t1(2);

tsp1(tsp1>t1) = [];
tsp2(tsp2>t1) = [];

tsp1 = tsp1-t0;
tsp2 = tsp2-t0;
tsp1(tsp1<0) = [];
tsp2(tsp2<0) = [];

if isempty(tsp1) || isempty(tsp2)
    N12 = 0;
    tsp12 = [];
    return; 
end
maxT = max(max(tsp1), max(tsp2));
T = ceil(maxT/deltat)+1;
spk = false(2,T);

ind1 = ceil(tsp1/deltat)+1;
ind2 = ceil(tsp2/deltat)+1;

ind1(isnan(ind1)) = [];
ind2(isnan(ind2)) = [];

spk(1,ind1) = 1;
spk(2,ind2) = 1;

spk12 = spk(1,:).*spk(2,:);
N12 = sum(spk12);
temp = find(spk12);
tsp12 = temp*deltat;


