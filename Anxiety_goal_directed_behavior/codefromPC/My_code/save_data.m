function save_data(data, fileName)
%% save data

numSpk = sum(sum(data.binMat));
numTrial = length(data.tsp);
numBoots = size(data.psth.ci,2)-1; 

%% tsp
tsp = zeros(numSpk, 1);
trial_id = zeros(numSpk, 1);
k = 0;
for m=1:numTrial
    tmptsp = data.tsp{m};
    n = length(tmptsp);
    tsp((k+1):(k+n)) = tmptsp;
    trial_id((k+1):(k+n)) = m;
    k = k+n;
end

%% psth
Fs = data.Fs;
psthSim = data.sim.psth;
psthFit = data.psth.mod;
psthCI = reshape(data.psth.ci(:,1:end-1), [], 1); 
T = length(psthFit);
t = (1:T)'/Fs; 
tCI = repmat(t, numBoots, 1); 

%% self-recovery function 
srfSim = data.sim.srf; 
ttStar = (1:length(srfSim))'; 
srfFit = data.srf.mod(1:length(srfSim)); 
srfCI = reshape(data.srf.ci(1:length(srfSim),1:end-1), [], 1); 
ttStarCI = repmat(ttStar, numBoots, 1); 

%% oscillatory 
phi = data.osc.hist.x'; 
oscSim = data.osc.hist.true'; 
oscFit = data.osc.mod; 
oscCI = reshape(data.osc.ci(:,1:end-1), [],1); 
phiCI = repmat(reshape(phi, [], 1), numBoots,1); 

%% save
save(fileName, 'tsp', 'trial_id', 't', 'Fs', 'psthSim', ...
    'psthFit', 'psthCI', 'tCI', 'ttStar', 'ttStarCI', ...
    'srfSim', 'srfFit', 'srfCI', 'phi', 'oscSim', 'oscFit', 'oscCI', 'phiCI');
