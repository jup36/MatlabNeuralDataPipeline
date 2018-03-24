function data = ci(data, Nboots)
%% generate confidence interval of fitted data with bootsmap method
if ~exist('data', 'var')
    disp('input arguments are not enought!');
    return;
end

if ~exist('Nboots', 'var')||isempty(Nboots)
    Nboots = 100; %number of bootstrap fitting
end

%% global parameters
numTrial = length(data.tsp);
% t0t1 = data.t0t1;
% t0 = t0t1(1);
% t1 = t0t1(2);
Fs = data.Fs;
% T = ceil((t1-t0)*Fs);

%% initial parameters for fitted variables
fitVar = data.fitVar;
for m=1:length(fitVar)
    if strcmp(fitVar{m}, 'psth')
        psthSim = data.psth.mod*Fs;
        data.psth.ci = zeros(length(psthSim), Nboots+1);
        data.psth.ci(:,end) = data.psth.mod;
    elseif strcmp(fitVar{m}, 'srf')
        lambda2 = data.srf.mod;
        data.srf.ci = zeros(length(lambda2), Nboots+1);
        data.srf.ci(:,end) = lambda2;
    elseif strcmp(fitVar{m}, 'osc')
        oscX = data.osc.basis.phi;
        oscY = data.osc.mod;
        oscFun = @(x, oscX, oscY)(interp1(oscX, oscY, x));
        oscPhase = reshape(data.oscPhase, [], 1);
        lambda3 = reshape(oscFun(oscPhase, oscX, oscY), [], numTrial);
        data.osc.ci = zeros(length(oscY), Nboots+1);
        data.osc.ci(:,end) = oscY;
    end
end
if ~exist('psthSim', 'var')
    try
        psthSim = data.sim.psthSim;
    catch
        psthSim = data.psth.mod*Fs;
    end
end
if ~exist('lambda2', 'var')
    try
        lambda2 = data.sim.lambda2;
    catch
        try
            lambda2 = data.srf.mod;
        catch
            lambda2 = 1;
        end
    end
end
if ~exist('oscFun', 'var')
    lambda3 = ones(1, numTrial);
end

%% start simulation
% tsp = cell(numTrial*Nboots,1);
% k = 0;
% for m=1:Nboots
%     for n=1:numTrial
%         k = k+1;
%         tsp{k} = simIMI(psthSim.*lambda3(:, n), lambda2);
%     end
% end

%% begin bootstrap
tic;
for m=1:Nboots
    for n=1:numTrial
        data.tsp{n} = simIMI(psthSim.*lambda3(:, n), lambda2);
    end

    data = anaData(data);
    for n=1:length(fitVar)
        if strcmp(fitVar{n}, 'psth')
            data.psth.ci(:,m) = data.psth.mod;
        elseif strcmp(fitVar{n}, 'srf')
            data.srf.ci(:,m) = data.srf.mod;
        elseif strcmp(fitVar{n}, 'osc')
            data.osc.ci(:,m) = data.osc.mod;
        end
    end
    if mod(m,10)==0
        toc;
        disp(m);
        tic;
    end
end
toc;




















