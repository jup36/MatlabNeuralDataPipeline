%% initialize parameters
clear data;
close all; 
clc; 

pars.tau = 10;      %time constant of the membrane 
pars.sig = 3*sqrt(pars.tau*2);       %amplitude of noise 
pars.vth = 20;      %voltage threshold
pars.vre = 0;       %resetting voltage
pars.i0 = 16;       %direct current to neuron 
pars.dt = 0.1;      %time step for simulation 
pars.total = 3000;  %total simulation time 
pars.taur = 1;      %refractory period

Fs = 1000;          %sampling rate 
dt = pars.dt;       %time step for simulation 
freq = 3;           %slow varying rate for input signal 
total = pars.total; %total simulation time 
data.t0t1 = [0,total/Fs];   %start and ending time point
data.pars = pars;   %parameters for simulation 
t = dt:dt:total;     %time 
% data.t = t;

amp = 1;            %amplitude of slow varying input signal 
temp = gen_ar(total*2, freq);   %generate AR signal at given frequency 
temp = (temp-mean(temp))/std(temp); %normalize and center
iin = amp*resample(temp, ceil(1/dt), 1);    %resample at simulation resolution 
tmp_n1 =ceil((length(iin)-length(t))/2);    
iin = iin(tmp_n1:(tmp_n1+length(t)-1)); %use data in center part to remove boundary effect

data.dc = iin+pars.i0;  %input signal for LIF neuron, this part is the same over all trials
% 
% figure;
% subplot(3,1,1:2);
% hold on;
N = 120;
data.tsp = cell(N,1);

sample_N = floor(length(t)*dt);  %since lfp is at resolution of dt, we want to down-sample it 
if sample_N/dt>length(t)
    sample_flag = 1;
end

data.lfp = zeros(sample_N, N);
lfp_freq = 40;          %each trial has its own signal at special frequency 
pars.lfp_freq = lfp_freq; 
lfp_amp =  1*ones(N,1); %linspace(0,2,N); %
pars.lfp_amp = lfp_amp; 
phase0 = rand(N,1); 
t1 = (1:(total*2))'/Fs; 

for m=1:N
    temp = gen_ar(total*2, lfp_freq);   %generate AR signal at special frequency
    temp = (temp-mean(temp))/std(temp);
    temp = FilterLFP([t1,temp], 'passband', [lfp_freq-5, lfp_freq+5]); 
    
    lfp = resample(temp(:,2), ceil(1/dt), 1);
    tmp_n1 = ceil((length(lfp)-length(t))/2);
    lfp = lfp_amp(m)*lfp(tmp_n1:(tmp_n1+length(t)-1)); %extract the centering part only
    
%     lfp = lfp_amp(m)*cos(2*pi*(t'*lfp_freq/1000+phase0(m))); 
    tsp = simLIF(pars,iin+lfp);     %simulation 
    data.tsp{m} = tsp;              %save spike trains
%     plot(tsp, (tsp>0)*m, '.');      %raster plot
    
    %save lfp
%     if sample_flag==1               %down-sampling
%         lfp(end+1:sample_N/dt) = lfp(end);
%     end
    temp = mean(reshape(lfp, 1/dt,[]));
    
    data.lfp(:,m) = temp;           %save lfp data 
end
data.trial_id = 1:N;
% xlim([0, total/Fs]);
% 
% subplot(313);
% plot(t/Fs, iin);

data.oscPhase = angle(hilbert(data.lfp)); 
data.Fs = 1000; 
data.maxISI = 50; 

data = anaData(data); 
% data.fitVar = {'psth', 'osc'}; 
    
    
    
    
    
    
    
    
    
    
    
    
    
    


