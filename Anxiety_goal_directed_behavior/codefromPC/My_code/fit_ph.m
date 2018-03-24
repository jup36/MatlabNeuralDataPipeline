function data = fit_ph(data)
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

%historyx
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

%binary vector for spike train
if ~isfield(data, 'binMat')
    data.binMat = tsp2bin(data.tsp, t0t1(1), t0t1(2), data.Fs);
end
Y = reshape(full(data.binMat), [], 1);

%%
%fit psth+oscillation first
X = X1+(rand(size(X1))-0.5)*0.00001;
% try
%     b_p = irls(X, Y, data.psth.beta);
% catch
    b_p = irls(X, Y);
% end
data.psth.beta = b_p;
data.psth.mod = exp(xbeta(data.psth.basis.basisX, data.psth.beta));
offset_p = xbeta(X1, b_p);
data.offset_p = offset_p;

%fit history
X = X2+(rand(size(X2))-0.5)*0.00001;
% try
%     b_h = irls(X, Y, data.srf.beta, offset_p);
% catch
    b_h = irls(X, Y, [], offset_p);
% end
data.srf.beta = b_h;
data.srf.mod = exp(xbeta(data.srf.basis.basisX, data.srf.beta));
data.offset = offset_p+xbeta(X2, b_h);
