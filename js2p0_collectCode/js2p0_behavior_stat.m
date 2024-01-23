
filePath = {'/Volumes/dudmanlab/junchol/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22)
            '/Volumes/dudmanlab/junchol/js2p0/WR37_022219', ... % Dual recording without silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked       
            '/Volumes/dudmanlab/junchol/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked      
            '/Volumes/dudmanlab/junchol/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/dudmanlab/junchol/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked

figSaveDir = '/Volumes/dudmanlab/junchol/js2p0/collectData/collectFigure/block_behavior'; 
        
%%
success_rate = cell(length(filePath),8); 
rch_angle = cell(length(filePath),8); 
pull_force = cell(length(filePath),8); 
rch_tort = cell(length(filePath),8); 

for f = 1:length(filePath)
    cd(filePath{f})
    file = dir('**/*js2p0_tbytSpkHandJsTrjBin_WR*.mat');
    load(fullfile(file.folder, file.name),'ss')
    b_id = [ss.blNumber]; 
    
    trType = {ss.trialType}; % trial type 
    rchDeg = {ss.rchAngDeg}; % reach angle in degrees
    maxForce = {ss.maxPullForce}; % max pull force  
    tort = {ss.hXYtort}; % tortuosity 
    
    for b = unique(b_id)
        b_idx = [ss.blNumber]==b;    
        % success rate
        success_rate{f,b} = cell2mat(cellfun(@(a) strcmpi(a, 'sp'), trType(b_idx), 'un', 0))';
        oms_rate{f,b} = cell2mat(cellfun(@(a) strcmpi(a, 'to'), trType(b_idx), 'un', 0))';
        success_rate_wo_oms{f,b} = success_rate{f,b}(~oms_rate{f,b}); 
        
        if b == 1 
           success_trials = find(success_rate{f,b}); 
           success_rate{f,b} = success_rate{f,b}(success_trials(3):end); 
        end 
        % reach angle
        rch_angle{f,b} = cell2mat(rchDeg(b_idx))';
        rch_angle{f,b} = rch_angle{f,b}(rch_angle{f,b} ~= 0); % get rid of zeros
        % max pull force
        pull_force{f,b} = cell2mat(maxForce(b_idx))'; 
        pull_force{f,b} = pull_force{f,b}(pull_force{f,b} ~= 0); % get rid of zeros 
        % tortuosity
        rch_tort{f,b} = cell2mat(tort(b_idx))';
        rch_tort{f,b} = rch_tort{f,b}(rch_tort{f,b} ~= 15); % 15 is the saturation value       
    end
    fprintf('processed file #%d\n', f) 
end

%save(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'))
load(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'))

%% descriptive stats and blockwise plot
% reach angle
avg_rch_angle = cell2mat(cellfun(@nanmean, rch_angle, 'un', 0)); 
blockwiseAvgPlot(avg_rch_angle)
%print(fullfile(figSaveDir,'blockwise_rchAngle'), '-dpdf', '-painters')

% reach angle median subtraction
rch_angle_med_sub = {}; 
for j = 1:size(rch_angle, 1)
    med_rch_angle = nanmedian(cell2mat(rch_angle(j, :)')); 
    rch_angle_med_sub = [rch_angle_med_sub; cellfun(@(a) a-med_rch_angle, rch_angle(j, :), 'un', 0)]; 
end
avg_rch_angle_med_sub = cell2mat(cellfun(@nanmean, rch_angle_med_sub, 'un', 0)); 
blockwiseAvgPlot(avg_rch_angle_med_sub)
ylim([-50 40])
%print(fullfile(figSaveDir,'blockwise_rchAngle_med_sub'), '-dpdf', '-painters')

% reach tortuosity
avg_rch_tort = cell2mat(cellfun(@nanmean, rch_tort, 'un', 0)); 
blockwiseAvgPlot(avg_rch_tort)
ylim([0 7])
print(fullfile(figSaveDir,'blockwise_rchTort'), '-dpdf', '-painters')
% pull force
avg_pull_force = -cell2mat(cellfun(@nanmean, pull_force, 'un', 0)); 
blockwiseAvgPlot(avg_pull_force)
print(fullfile(figSaveDir,'blockwise_pullForce'), '-dpdf', '-painters')
% success rate
avg_s_rate = cell2mat(cellfun(@(a) sum(a)/length(a), success_rate, 'un', 0)).*100; 
blockwiseAvgPlot(avg_s_rate)
print(fullfile(figSaveDir,'blockwise_successRate'), '-dpdf', '-painters')
% success rate excluding omissions
avg_s_rate_wo_oms = cell2mat(cellfun(@(a) sum(a)/length(a), success_rate_wo_oms, 'un', 0)).*100; 
blockwiseAvgPlot(avg_s_rate_wo_oms)
ylim([20 110])
print(fullfile(figSaveDir,'blockwise_successRate_without_omission'), '-dpdf', '-painters')
% omission rate
avg_oms_rate = cell2mat(cellfun(@(a) sum(a)/length(a), oms_rate, 'un', 0)).*100; 
blockwiseAvgPlot(avg_oms_rate)
print(fullfile(figSaveDir,'blockwise_omissionRate'), '-dpdf', '-painters')


%% plot 
function blockwiseAvgPlot(m_b_mat)
    cmat = {'b', 'b', 'r', 'r', 'b', 'b', 'r', 'r', 'b', 'b', 'r', 'r', 'b', 'b', 'r', 'r'}; 
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
                scatter(x{m,b}, y{m,b}, 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', cmat{b}, ...
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
    ypad = floor((max(m_b_mat(:))-min(m_b_mat(:)))*0.08);
    ylim([min(m_b_mat(:))-ypad, max(m_b_mat(:))+ypad])
end



