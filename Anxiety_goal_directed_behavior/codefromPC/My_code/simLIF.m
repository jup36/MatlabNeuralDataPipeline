function [tsp, v, pars] = simLIF(pars, iin)
%% simulate LIF neuron: tau*dV/dt = -V+I, V(t)=Vth-> V(t+dt) = V0
%input:
%   iin: input current
%   pars: structure data containing useful parameters
%       %sig: effect of white noise, default 0
%       %tau: time constant of neuron, default 1
%       %vth: threshold for firing spike, default 0.5
%       %vre: reseting threshold, default 0
%       %taur: refractory period, default 0
%       %dt:  time step for simulation, default 0.001
%       %total: total simulation time
%output:
%   tsp: time of stamp
%   v: voltage trace
%   t: time

%% initialize
Fs = 1000; 
if ~exist('pars', 'var')
    i0 = 0.5;
    sig = 0;
    tau = 1;
    vth = 0.5;
    vre = 0;
    taur = 0;
    dt = 0.001;
    total = 1000;
else
    i0 = pars.i0;
    try sig = pars.sig;
    catch
        sig = 0; 
    end
    tau = pars.tau;
    vth = pars.vth;
    vre = pars.vre;
    try taur = pars.taur;
    catch
        taur = 0;
    end
    dt = pars.dt;
    total = pars.total;
end
Ntotal = ceil(total/dt);

%total simulation time
if ~exist('iin','var')
    iin = zeros(Ntotal,1);
elseif  length(iin)==1
    iin = ones(Ntotal,1)*iin;
elseif length(iin)>=Ntotal
    iin(Ntotal+1:end) = [];
else
    Ntotal = length(iin);
    pars.total = Ntotal*dt;
end
iin = iin+i0;

v = ones(size(iin))*vre;
v(1) = vre+exp(rand(1)*log(vth-vre)); 
w = randn(size(iin));

%% start simulation
count_refrac = 0;
count_spk = 0;
tsp = zeros(ceil(Ntotal*dt/tau),1);

dttau = dt/tau;
sqrtdttau = sqrt(dt)/tau;

for m=2:Ntotal
    if count_refrac==0
        v(m) = v(m-1)*(1-dttau)+dttau*iin(m)+sqrtdttau*w(m)*sig;
    else
        count_refrac = count_refrac-1;
        continue;
    end
    
    if v(m)>vth
        count_refrac = ceil(taur/dt);
        v(m) = vre;
        count_spk = count_spk+1;
        tsp(count_spk) = m*dt/Fs;
    end
end

tsp(count_spk+1:end) = [];


























