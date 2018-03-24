function data = anaData(data, hist_flag)
%% analyze data
%hist_flag: report histogram method or not

if ~exist('data', 'var')
    data = evalin('base', 'data');
end
if ~exist('hist_flag', 'var')
    hist_flag = 0; 
end
if ~isfield(data, 'fitVar')
    data.fitVar = {'psth', 'srf', 'osc'};
end
%% pre-processing
temp = data.t0t1;
t0 = temp(1);
t1 = temp(2);
Fs = data.Fs;
T = ceil((t1-t0)*Fs);
numTrial = length(data.tsp);

data.binMat = tsp2bin(data.tsp, t0, t1, Fs);

%spline basis for psth
if any(strcmp(data.fitVar, 'psth'))
    if ~isfield(data, 'psth')
        if ~isfield(data, 'psth_flag')
            data.psth_flag = 0; 
        end
        data.psth.basis = nsPsth(data, data.psth_flag);
        data.psth.type = 0;
        data.psth.xind = repmat((1:T)', numTrial, 1);
    end
    
    if isfield(data.psth, 'X')
        data.psth = rmfield(data.psth, 'X');
    end
end
%spline basis for self-recovery function
if any(strcmp(data.fitVar, 'srf'))
    try
        maxISI = data.maxISI;
    catch
        maxISI = 100;
        data.maxISI = maxISI; 
    end

    if ~isfield(data, 'srf')
        if ~isfield(data, 'srf_flag')
            data.srf_flag = 0; 
        end
        data.srf.basis = nsSrf(data, data.srf_flag);
        data.srf.type = 1;
    end
    
    tAfter = ttStar(data.tsp, t0, t1, 1, Fs);
    xind = reshape(tAfter, [], 1);
    data.srf.xind = xind;
    data.srf.ind = ((xind>=1)&(xind<=round(maxISI*Fs/1000)));  
end
%spline basis for oscillation
if any(strcmp(data.fitVar, 'osc'))
    if ~isfield(data, 'osc')
        if ~isfield(data, 'osc_flag')
            data.osc_flag = 0; 
        end
        data.osc.basis = nsOsc(data, data.osc_flag);
        data.osc.type = 1;
    end
    
    dphi = diff(data.osc.basis.phi(1:2));
    phi0 = min(data.osc.basis.phi);
    data.osc.xind = max(ceil(reshape(data.oscPhase-phi0,[], 1)/dphi), 1);
end

% spline basis for network effect 
if any(strcmp(data.fitVar, 'net'))
    try 
        maxN = data.maxN; 
    catch 
        maxN = 1; 
        data.maxN = maxN; 
    end
    
    if ~isfield(data, 'net')
        if ~isfield(data, 'net_flag')
            data.net_flag = 0; 
        end
        data.net.basis = nsNet(data, data.net_flag);
        data.net.type = 1;
    end
    
    dn = diff(data.net.basis.n(1:2));
    n0 = min(data.net.basis.n);
    temp = max(ceil(reshape(data.netN-n0,[], 1)/dn), 1); 
    temp(data.netN(:)>maxN) = length(data.net.basis.n); 
    data.net.xind = temp;
    data.net.ind = true(size(temp)); 
end
% 
% if any(strcmp(data.fitVar, 'net'))
%     if ~isfield(data, 'net')
%         if ~isfield(data, 'net_flag')
%             data.net_flag = 0; 
%         end
%         data.net.basis = nsNet(data, data.net_flag); 
%         data.net.type = 1; 
%         data.net.xind = reshape(data.sum_spk, [], 1); 
%         data.net.ind = (data.net.xind>=1) & (data.net.xind<=data.net.basis.n(end)); 
%     end
% end
% if any(strcmp(data.fitVar, 'net'))
%     try 
%         maxN = data.maxN; 
%     catch 
%         maxN = 1; 
%         data.maxN = maxN; 
%     end
%     
%     if ~isfield(data, 'net')
%         if ~isfield(data, 'net_flag')
%             data.net_flag = 0; 
%         end
%         data.net.basis = nsNet(data, data.net_flag);
%         data.net.type = 1;
%     end
%     
%     dn = diff(data.net.basis.n(1:2));
%     n0 = min(data.net.basis.n);
%     temp = max(ceil(reshape(data.sum_spk-n0,[], 1)/dn), 1); 
%     temp(temp>maxN) = length(data.net.basis.n); 
%     data.net.xind = temp;
%     data.net.ind = (data.srf.xind<10); 
% end
%% fit GLM model

% if isfield(data, 'offset')
%     data = fitData(data, data.offset);
% else
if isfield(data, 'penalty')
    data = fitData(data, data.penalty);
else
    data = fitData(data);
end
% end


%% histogram results
if hist_flag==1
    if any(strcmp(data.fitVar, 'psth'))
        data.psth.hist = getHistPsth(data.tsp, t0, t1);
    end
    if any(strcmp(data.fitVar, 'srf'))
        data.srf.hist = getHistSrf(data.tsp, t0, t1, Fs, maxISI);
    end
    if any(strcmp(data.fitVar, 'osc'))
        data.osc.hist = getHistOsc(data.oscPhase, data.binMat, linspace(-pi, pi, 100));
    end
end
