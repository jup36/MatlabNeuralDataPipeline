function c12 = sfc(tsp, lfp, params, t0t1)
%% compute spike field coherence

if nargin<2
    disp('input arguments are not enough');
end

%params
if ~exist('params', 'var')||isempty('params')
    params = struct();
end
if ~isfield(params, 'Fs')
    params.Fs = 1000;
end
if ~isfield(params, 'fpass')
    params.fpass = [0, 100];
end
if ~isfield(params, 'tapers')
    params.tapers = [5, 9];
end
if ~isfield(params, 'err')
    params.err = [1, 0.05];
end
if ~isfield(params, 'trialave')
    params.trialave = 1; 
end 
Fs = params.Fs;

%time period
if exist('t0t1', 'var')
    ind0 = max(ceil(t0t1(1)*Fs), 1);
    ind1 = max(ceil(t0t1(2)*Fs), 2);
    lfp = double(lfp(ind0:ind1, :));
else
    ind0 = 1;
    ind1 = size(lfp,2);
    t0t1 = [ind0-1, ind1]/Fs;
end

%use chronux to calculate 
datasp = format_tsp(tsp, t0t1);
datalfp = locdetrend(lfp, Fs, [0.1, 0.05]); 
[C, phi, S12, S1, S2, f, zerosp, confC, phistd] = coherencycpt(datalfp, datasp, params);
c12.C = C;
c12.phi = phi;
c12.S12 = S12;
c12.S1 = S1;
c12.S2 = S2;
c12.f = f;
c12.zerosp = zerosp;
c12.confC = confC;
c12.phistd = phistd;

function datasp = format_tsp(tsp, t0t1)
%% change tsp to default format of chronux
numTrial = length(tsp);
datasp(numTrial) = struct;
if exist('t0t1', 'var')
    t0 = t0t1(1);
    t1 = t0t1(2);
else
    t0 = -1;
    t1 = inf;
end

for m=1:numTrial
    temp = tsp{m};
    temp(temp<=t0) = [];
    temp(temp>t1) = [];
    datasp(m).time = temp-t0;
end

















