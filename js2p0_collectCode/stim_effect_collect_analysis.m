
filePath_cg = {'/Volumes/Extreme SSD/js2p0/WR37_022619/Matfiles', ...  % Cg recording contra-Cg silencing, Trj checked
               '/Volumes/Extreme SSD/js2p0/WR38_052319/Matfiles', ...  % Cg recording contra-Cg silencing, Trj checked
               '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles'};     % Dual recording with contra Cg silencing, Trj checked
                
filePath_delayed_cg = {'/Volumes/Extreme SSD/js2p0/WR37_022719/Matfiles', ... % B Only contra-Cg delayed silencing, Trj checked
                       '/Volumes/Extreme SSD/js2p0/WR39_091019/Matfiles', ... % B Only contra-Cg delayed silencing, Trj checked 
                       '/Volumes/Extreme SSD/js2p0/WR39_091119/Matfiles', ... % B Only contra-Cg delayed silencing, Trj checked 
                       '/Volumes/Extreme SSD/js2p0/WR39_100319/Matfiles', ... % B Only contra-Cg delayed silencing, Trj checked      
                       '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles'};    % Dual recording with contra Cg delayed silencing, Trj checked 

filePath_m1 = {'/Volumes/Extreme SSD/js2p0/WR38_050119/Matfiles', ...  % M1 silencing
               '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ...  % Corticostriatal recording M1 silencing, Trj checked     
               '/Volumes/Extreme SSD/js2p0/WR40_082219/Matfiles', ...
               '/Volumes/Extreme SSD/js2p0/WR45_030220', ...
               '/Volumes/Extreme SSD/js2p0/WR46_012520'};     % M1 silencing

%% Cg stim without delay
reach_prob.cg_no_delay.dat = nan(length(filePath_cg), 2); 

for j = 1:length(filePath_cg)
    file_dir = dir(fullfile(filePath_cg{j}, 'stim_effect_rez_*'));  
    load(fullfile(file_dir.folder, file_dir.name), 'jsReachOn_prob_noStim', 'jsReachOn_prob_stim');
    reach_prob.cg_no_delay.dat(j, 1) = jsReachOn_prob_noStim; 
    reach_prob.cg_no_delay.dat(j, 2) = jsReachOn_prob_stim; 
    clearvars jsReachOn_prob_noStim jsReachOn_prob_stim
end

% paired t-test
[~, reach_prob.cg_no_delay.p, ~, reach_prob.cg_no_delay.stats] = ttest(reach_prob.cg_no_delay.dat(:, 1), reach_prob.cg_no_delay.dat(:, 2)); 
reach_prob_scatter_plot(reach_prob.cg_no_delay.dat)
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'reach_prob_cg_no_delay'), '-painters', '-dpdf')
    
%% Cg delayed stim
reach_prob.cg_with_delay.dat = nan(length(filePath_delayed_cg), 2); 

for j = 1:length(filePath_delayed_cg)
    file_dir = dir(fullfile(filePath_delayed_cg{j}, 'stim_effect_rez_*'));  
    load(fullfile(file_dir.folder, file_dir.name), 'jsReachOn_prob_noStim', 'jsReachOn_prob_stim');
    reach_prob.cg_with_delay.dat(j, 1) = jsReachOn_prob_noStim; 
    reach_prob.cg_with_delay.dat(j, 2) = jsReachOn_prob_stim; 
    clearvars jsReachOn_prob_noStim jsReachOn_prob_stim
end

% paired t-test
[~, reach_prob.cg_with_delay.p, ~, reach_prob.cg_with_delay.stats] = ttest(reach_prob.cg_with_delay.dat(:, 1), reach_prob.cg_with_delay.dat(:, 2)); 
reach_prob_scatter_plot(reach_prob.cg_with_delay.dat)
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'reach_prob_cg_with_delay'), '-painters', '-dpdf')

[~, reach_prob.cg_with_or_without_delay.p, ~, reach_prob.cg_with_or_without_delay.stats] = ...
    ttest([reach_prob.cg_no_delay.dat(:, 1); reach_prob.cg_with_delay.dat(:, 1)], [reach_prob.cg_no_delay.dat(:, 2); reach_prob.cg_with_delay.dat(:, 2)]); 

% Wilcoxon signed-rank test (for paired data)
[reach_prob.cg_with_or_without_delay.p_signedrank, ~, reach_prob.cg_with_or_without_delay.stats_signedrank] = signrank([reach_prob.cg_no_delay.dat(:, 1); reach_prob.cg_with_delay.dat(:, 1)], [reach_prob.cg_no_delay.dat(:, 2); reach_prob.cg_with_delay.dat(:, 2)], 'tail', 'right');

reach_prob_scatter_plot([reach_prob.cg_no_delay.dat; reach_prob.cg_with_delay.dat])
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'reach_prob_cg_with_and_without_delay'), '-painters', '-dpdf')
    
%% M1 stim
reach_prob.m1.dat = nan(length(filePath_m1), 2); 

for j = 1:length(filePath_m1)
    file_dir = dir(fullfile(filePath_m1{j}, 'stim_effect_rez_*'));  
    load(fullfile(file_dir.folder, file_dir.name), 'jsReachOn_prob_noStim', 'jsReachOn_prob_stim');
    reach_prob.m1.dat(j, 1) = jsReachOn_prob_noStim; 
    reach_prob.m1.dat(j, 2) = jsReachOn_prob_stim; 
    clearvars jsReachOn_prob_noStim jsReachOn_prob_stim
end

% paired t-test
[~, reach_prob.m1.p, ~, reach_prob.m1.stats] = ttest(reach_prob.m1.dat(:, 1), reach_prob.m1.dat(:, 2)); 
% Wilcoxon signed-rank test (for paired data)
[reach_prob.m1.p_signedrank, ~, reach_prob.m1.stats_signedrank] = signrank(reach_prob.m1.dat(:, 1), reach_prob.m1.dat(:, 2), 'tail', 'right');
% Wilcoxon rank-sum test (for independent data)
%[reach_prob.m1.p_ranksum, ~, reach_prob.m1.stats_ranksum] = ranksum(reach_prob.m1.dat(:, 1), reach_prob.m1.dat(:, 2));

reach_prob_scatter_plot(reach_prob.m1.dat)
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'reach_prob_m1_updated'), '-painters', '-dpdf')
    
%% Independent t-test between M1 vs. Cg
reach_prob.cg_no_delay.stim_effect_size = reach_prob.cg_no_delay.dat(:, 2)-reach_prob.cg_no_delay.dat(:, 1);  
reach_prob.cg_with_delay.stim_effect_size = reach_prob.cg_with_delay.dat(:, 2)-reach_prob.cg_with_delay.dat(:, 1);  
reach_prob.m1.stim_effect_size = reach_prob.m1.dat(:, 2)-reach_prob.m1.dat(:, 1); 

[~, reach_prob.cg_vs_m1.p, ~, reach_prob.cg_vs_m1.stats] = ttest2([reach_prob.cg_no_delay.stim_effect_size; reach_prob.cg_with_delay.stim_effect_size], reach_prob.m1.stim_effect_size); 

%save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'reach_prob_rez_stats.mat'),'reach_prob')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reach_prob_scatter_plot(m_b_mat)
    cmat = {'k', 'b'}; 
    randX = [-.5 + rand(size(m_b_mat,1)*size(m_b_mat,2),1)].*0.25; 
    alph = [0.8, 0.4]; 
    x = cell(size(m_b_mat,1), size(m_b_mat,2)); 
    y = cell(size(m_b_mat,1), size(m_b_mat,2)); 
    
    figure; hold on
    val = 1; 
    for b = 1:size(m_b_mat,2)
        for m = 1:size(m_b_mat,1)
            x{m,b} = b+randX(val); 
            y{m,b} = m_b_mat(m,b); 
            if ~isnan(y{m,b})
                scatter(x{m,b}, y{m,b}, 100, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmat{b}, ...
                'MarkerFaceAlpha', alph(mod(b,2)+1))
                % connect dots across blocks
                if b > 1
                   plot([x{m,b-1}, x{m,b}], [y{m,b-1}, y{m,b}], 'k:')    
                end
                val = val + 1; 
            end
        end
        avg_y_b = nanmean([y{:,b}]); 
        % mark the block average
        plot([b-.3, b+.3], [avg_y_b, avg_y_b], cmat{b}, 'LineWidth', 2) 
    end
    hold off
    set(gca, 'TickDir', 'out')
    
    xlim([.5 size(m_b_mat,2)+.5])
    ylim([0 .7])
end