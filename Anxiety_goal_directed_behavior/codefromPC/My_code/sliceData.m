function newdata = sliceData(data, ind, t0t1)
%% extract specific trials or time interval
%ind: trials to be selected
%t0t1: time interval

tsp = data.tsp;
numTrial = length(tsp);
useTrial = false(numTrial, 1);

if ~exist('ind', 'var')||isempty(ind)
    useTrial = true(size(useTrial));
elseif islogical(ind)
    useTrial = ind;
else
    ind = unique(ind);
    ind(ind>numTrial) = [];
    ind(ind<1) = [];
    useTrial(ind) = 1;
end

newdata.Fs = data.Fs;
newdata.tsp = tsp(useTrial);

if exist('t0t1', 'var')
    newdata.t0t1 = t0t1-t0t1(1);
    for m=1:length(newdata.tsp)
        temp = newdata.tsp{m};
        temp(temp<=t0t1(1)) = [];
        temp(temp>t0t1(2)) = [];
        newdata.tsp{m} = temp-t0t1(1);
    end
else
    newdata.t0t1 = data.t0t1;
    t0t1 = data.t0t1;
end
ind0 = floor(t0t1(1)*data.Fs)+1;
ind1 = ceil(t0t1(2)*data.Fs);
newdata.maxISI = data.maxISI;

if isfield(data, 'sim')
    newdata.sim = data.sim;
end

if isfield(data, 'lfp')
    newdata.lfp = data.lfp(ind0:ind1, useTrial);
end

if isfield(data, 'oscPhase')
    newdata.oscPhase = data.oscPhase(ind0:ind1, useTrial);
end

if isfield(data, 'binMat')
    newdata.binMat = data.binMat(ind0:ind1, useTrial);
end

if isfield(data, 'fitVar')
    newdata.fitVar = data.fitVar;
end

if isfield(data, 'psth')
    temp = data.psth;
    if isfield(temp, 'xind')
        newdata.psth = rmfield(newdata.psth, 'xind');
    end
    
    if isfield(temp, 'X')
        newdata.psth = rmfield(newdata.psth, 'X');
    end
end

if isfield(data, 'srf')
    newdata.srf = data.srf;
    if isfield(newdata.srf, 'xind')
        newdata.srf = rmfield(newdata.srf, 'xind');
    end
    
    if isfield(newdata.srf, 'X')
        newdata.srf = rmfield(newdata.srf, 'X');
    end
end

if isfield(data, 'osc')
    newdata.osc = data.osc;
    if isfield(newdata.osc, 'xind')
        newdata.psth = rmfield(temp, 'xind');
    end
    
    if isfield(newdata.osc, 'X')
        newdata.psth = rmfield(newdata.psth, 'X');
    end
end
