function fig = plotData(data, alpha)
%% plot fitted data
%alpha: (1-alpha)*100% confidence band
if ~exist('alpha', 'var')
    alpha = 0.05;
end

fig = figure;
num_trial = length(data.tsp);
t0t1 = data.t0t1;
Fs = data.Fs;
T = ceil(diff(t0t1)*Fs);
t = t0t1(1)+(1:T)'/Fs;

%% raster plot
subplot(221); hold on;
for m=1:num_trial
    tmp_tsp = data.tsp{m};
    plot(tmp_tsp, m*ones(size(tmp_tsp)), '.', 'markersize', 3);
end
xlim(t0t1);
ylim([0, num_trial]);
set(gca, 'xtick', []);
ylabel('Trial #', 'fontsize', 15);
title('Raster Plot', 'fontsize', 15);

%% psth
subplot(223); hold on;

if isfield(data.psth, 'ci')
    tmp_sd = std(data.psth.ci(:, 1:(end-1)), 0, 2);
    y1 = data.psth.ci(:, end)+tmp_sd*norminv(1-alpha/2);
    y2 = data.psth.ci(:, end)-tmp_sd*norminv(1-alpha/2);
    
    %     y1 = quantile(data.psth.ci(:, 1:(end-1)), alpha/2, 2);
    %     y2 = quantile(data.psth.ci(:, 1:(end-1)), 1-alpha/2, 2);
    y = [y1; flipud(y2)];
    t1 = [t+rand(size(t))/Fs/100; flipud(t)];
    fill(t1, y*Fs, 'c', 'edgecolor', 'c');
    plot(t, data.psth.ci(:, end)*1000, 'linewidth', 3);
elseif isfield(data.psth, 'mod')
    plot(t, data.psth.mod*Fs, 'linewidth', 3);
end
xlim(t0t1);
xlabel('Time (Sec.)', 'fontsize', 15);
ylabel('Firing rate (Hz)', 'fontsize', 15);
title('PSTH', 'fontsize', 15);

%% self recovery function
subplot(222); hold on;
t = (1:data.maxISI)';
if isfield(data, 'srf')
    if isfield(data.srf, 'ci')
        tmp_sd = std(data.srf.ci(:, 1:(end-1)), 0, 2);
        y1 = data.srf.ci(:, end)+tmp_sd*norminv(1-alpha/2);
        y2 = data.srf.ci(:, end)-tmp_sd*norminv(1-alpha/2);
        %         y1 = quantile(data.srf.ci(:, 1:(end-1)), alpha/2, 2);
        %         y2 = quantile(data.srf.ci(:, 1:(end-1)), 1-alpha/2, 2);
        y = [y1; flipud(y2)];
        t1 = [t+rand(size(t))*0.0001; flipud(t)];
        fill(t1, y, 'c', 'edgecolor', 'c');
        plot(t, data.srf.ci(:, end), 'linewidth', 3);
    elseif isfield(data.srf, 'mod')
        plot(t, data.srf.mod, 'linewidth', 3);
    end
end
xlim([0, data.maxISI]);
xlabel('t-t^*', 'fontsize', 15);

% ylabel('');
title('Self-Recovery Function', 'fontsize', 15);

%% phase modulation
subplot(224); hold on;
if isfield(data, 'osc')
    phi = reshape(data.osc.basis.phi, [], 1);
    if isfield(data.osc, 'ci')
        tmp_sd = std(data.osc.ci(:, 1:(end-1)), 0, 2);
        y1 = data.osc.ci(:, end)+tmp_sd*norminv(1-alpha/2);
        y2 = data.osc.ci(:, end)-tmp_sd*norminv(1-alpha/2);
        
        %         y1 = quantile(data.osc.ci(:, 1:(end-1)), alpha/2, 2);
        %         y2 = quantile(data.osc.ci(:, 1:(end-1)), 1-alpha/2, 2);
        y = [y1; flipud(y2)];
        phi1 = [phi+rand(size(phi))*(0.00001); flipud(phi)];
        fill(phi1, y, 'c', 'edgecolor', 'c');
        plot(phi, data.osc.ci(:, end), 'linewidth', 3);
    elseif  isfield(data.osc, 'mod')
        plot(phi, data.osc.mod, 'linewidth', 3);
    end
end

xlim([-pi, pi]);
xlabel('Phase', 'fontsize', 15);
% ylabel('Firing rate (Hz)');
title('Phase Modulation', 'fontsize', 15);






















