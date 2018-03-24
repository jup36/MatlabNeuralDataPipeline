function data = fit_pho(data)
%% fit glm model of spike trains. three factors are included: psth, history, oscillation
t0t1 = data.t0t1;
% T = ceil(data.Fs*diff(t0t1));
num_trial = length(data.tsp);

%psth
if ~isfield(data, 'psth_flag')
    data.psth_flag = 1;
end

if ~isfield(data, 'psth') || ~isfield(data.psth, 'basis')
    data.psth.basis = nsPsth(data, data.psth_flag);
    data.psth_flag = data.psth.basis.knots;
end
X1 = repmat(data.psth.basis.basisX, num_trial, 1);

%history
if ~isfield(data, 'srf_flag')
    data.srf_flag = 0;
end
if ~isfield(data, 'maxISI')
    data.maxISI = 100;
end
if ~isfield(data, 'srf') || ~isfield(data.srf, 'basis')
    data.srf.basis = nsSrf(data, data.srf_flag);
    data.srf_flag = data.srf.basis.knots;
end
tAfter = ttStar(data.tsp, t0t1(1), t0t1(2), 1, data.Fs);
xind = reshape(tAfter, [], 1);
max_ind = ceil(data.maxISI*data.Fs/1000);
xind(xind>max_ind) = max_ind;
xind(isnan(xind)) = max_ind;
data.srf.xind = xind;
X2 = data.srf.basis.basisX(xind, :);

%oscillation
if ~isfield(data, 'osc_flag')
    data.osc_flag = 0;
end
if ~isfield(data, 'osc') || ~isfield(data.osc, 'basis')
    data.osc.basis = nsOsc(data, data.osc_flag);
    data.osc_field = data.osc.basis.knots;
end
dphi = diff(data.osc.basis.phi(1:2));
phi0 = min(data.osc.basis.phi);
data.osc.xind = max(ceil(reshape(data.oscPhase-phi0,[], 1)/dphi), 1);
X3 = data.osc.basis.basisX(data.osc.xind, :);

%binary vector for spike train
if ~isfield(data, 'binMat')
    data.binMat = tsp2bin(data.tsp, t0t1(1), t0t1(2), data.Fs);
end
Y = reshape(full(data.binMat), [], 1);

%%
n1 = size(X1, 2);
% n2 = size(X2, 2);
% n3 = size(X3, 2);

%fit psth+oscillation first
X = [X1, X3];
X = X+(rand(size(X))-0.5)*0.00001;
% try
%     b_po = irls(X, Y, data.b_po);
% catch
    b_po = irls(X, Y);
% end
temp = xbeta(data.osc.basis.basisX, b_po((n1+2):end));
c_o = -log(mean(exp(temp)));   %constant shift for oscillation
data.osc.beta = [c_o; b_po((n1+2):end)];
data.osc.mod = exp(xbeta(data.osc.basis.basisX, data.osc.beta));
c_p = b_po(1)-c_o;
data.psth.beta =[c_p; b_po(2:(n1+1))];
data.psth.mod = exp(xbeta(data.psth.basis.basisX, data.psth.beta));
offset_po = xbeta([X1, X3], b_po);
data.offset_po = offset_po;
data.b_po = b_po;

%fit history
X = X2+(rand(size(X2))-0.5)*0.00001;
% try
%     b_h = irls(X, Y, data.srf.beta, offset_po);
% catch
    b_h = irls(X, Y, [], offset_po);
% end
data.srf.beta = b_h;
data.srf.mod = exp(xbeta(data.srf.basis.basisX, data.srf.beta));
data.offset = offset_po+xbeta(X2, b_h);
