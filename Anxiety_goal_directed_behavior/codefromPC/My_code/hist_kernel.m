function ihbasis  = hist_kernel(maxISI)
%% make nonlinearly stretched basis consisting of raised cosines:
%%this is used for generating kernel of history effect
%Input:
%       n: number of basis
%       mxt: maximum influential time of one spike
%       dt: simulation time
%       taur: refractory period
%       b(optional) : offset for nonlinear stretching of x axis:  y = log(x+b)
%           (larger b -> more nearly linear stretching)
%Output:
%       iht: time lattice on which basis is defined
%       ihbas = orthogonalized basis
%       ihbasis = original basis (non-orthogonal)
%
%This code is a revised version of Jonathan Pillow.
%Pengcheng Zhou, zhoupc1988@gmail.com

b = 5; 
taur = 2; 
dt = 1; 
nb = 5; 

%% construct raised cosine basis
ctrb = logspace(log10(taur+b), log10(maxISI+b), nb-1);     %center of basis, equal distance in log space between taur and mxt
db = range(log(ctrb))/(nb-2);       %spacing between raised cosine peaks.
mxt = ctrb(end)*exp(2*db)-b;       %maximum time bin
iht = ceil((1:dt:mxt)');
nt = length(iht);       %number of points in iht
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector
ihbasis = ff(repmat(log(iht+b), 1, nb-1), repmat(log(ctrb), nt, 1), db);

% create first basis vector as step-function for absolute refractory period
ind = find(iht<=taur);
ih0 = zeros(size(ihbasis,1),1);
ih0(ind) = 1;
ihbasis(ind,:) = 0;
ihbasis = [ih0,ihbasis];


