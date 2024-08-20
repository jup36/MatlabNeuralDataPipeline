
filePath = {'/Volumes/Extreme SSD/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22)
            '/Volumes/Extreme SSD/js2p0/WR37_022219', ... % Dual recording without silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked       
            '/Volumes/Extreme SSD/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked      
            '/Volumes/Extreme SSD/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked

figSaveDir = '/Volumes/Extreme SSD/js2p0/collectData/collectFigure/block_behavior'; 
        
%%
success_rate = cell(length(filePath),8); 
rch_angle = cell(length(filePath),8); 
pull_force = cell(length(filePath),8); 
rch_tort = cell(length(filePath),8); 
rch_dur = cell(length(filePath),8);
pull_dur = cell(length(filePath),8);  
reaction_time = cell(length(filePath), 8); 

for f = 1:length(filePath)
    cd(filePath{f})
    file = dir('**/*js2p0_tbytSpkHandJsTrjBin_WR*.mat');
    load(fullfile(file.folder, file.name), 'ss', 'jkvt')
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
        % reach, pull duration 
        [rch_dur{f,b}, pull_dur{f,b}] = bByb_reach_pull_duration(ss(b_idx)); 
        % reaction time
        rt_blk{f,b} = bByb_reaction_time(jkvt(b_idx), ss(b_idx)); 
        
        [max_rch_spd{f,b}, max_rch_spd_x{f,b}, max_rch_spd_y{f,b}, max_rch_spd_z{f,b}] = bByb_maxRchSpeed(ss(b_idx)); 
    end
    fprintf('processed file #%d\n', f) 
end

%save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'))
%load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'))

%% descriptive stats and blockwise plot
%% reach angle
avg_rch_angle = cell2mat(cellfun(@nanmean, rch_angle, 'un', 0)); 
blockwiseAvgPlot(avg_rch_angle)
%print(fullfile(figSaveDir,'blockwise_rchAngle'), '-dpdf', '-vector')

% reach angle median subtraction
rch_angle_med_sub = {}; 
for j = 1:size(rch_angle, 1)
    med_rch_angle = nanmedian(cell2mat(rch_angle(j, :)')); 
    rch_angle_med_sub = [rch_angle_med_sub; cellfun(@(a) a-med_rch_angle, rch_angle(j, :), 'un', 0)]; 
end
avg_rch_angle_med_sub = cell2mat(cellfun(@nanmean, rch_angle_med_sub, 'un', 0)); 
blockwiseAvgPlot(avg_rch_angle_med_sub)
ylim([-50 40])
%print(fullfile(figSaveDir,'blockwise_rchAngle_med_sub'), '-dpdf', '-vector')

% repeated measures analysis of variance
rm_rez.rch_angle = blockByblock_rm(avg_rch_angle_med_sub, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]}); 
rch_angle_m_b_avg_C = table2cell(rm_rez.rch_angle.rm_tab); 
rch_angle_m_b_avg = cell2mat(rch_angle_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(rch_angle_m_b_avg); 
ylim([-20 15])
%print(fullfile(figSaveDir,'blockwise_rchAngle_med_sub_m_b_combined'), '-dpdf', '-vector')

%% reach tortuosity
avg_rch_tort = cell2mat(cellfun(@nanmean, rch_tort, 'un', 0)); 
blockwiseAvgPlot(avg_rch_tort)
ylim([0 7])
%print(fullfile(figSaveDir,'blockwise_rchTort'), '-dpdf', '-vector')

% repeated measures analysis of variance
rm_rez.rch_tortuosity = blockByblock_rm(avg_rch_tort, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]}); 
rch_tort_m_b_avg_C = table2cell(rm_rez.rch_tortuosity.rm_tab); 
rch_tort_m_b_avg = cell2mat(rch_tort_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(rch_tort_m_b_avg); 
ylim([2 6])
%print(fullfile(figSaveDir,'blockwise_tortuosity_m_b_combined'), '-dpdf', '-vector')

%% pull force
avg_pull_force = -cell2mat(cellfun(@nanmean, pull_force, 'un', 0)); 
blockwiseAvgPlot(avg_pull_force)
%print(fullfile(figSaveDir,'blockwise_pullForce'), '-dpdf', '-vector')

% repeated measures analysis of variance
rm_rez.pull_force = blockByblock_rm(avg_pull_force, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
pull_force_m_b_avg_C = table2cell(rm_rez.pull_force.rm_tab); 
pull_force_m_b_avg = cell2mat(pull_force_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(pull_force_m_b_avg); 
ylim([0 140])
%print(fullfile(figSaveDir,'blockwise_pull_force_m_b_combined'), '-dpdf', '-vector')

%% success rate
avg_s_rate = cell2mat(cellfun(@(a) sum(a)/length(a), success_rate, 'un', 0)).*100; 
blockwiseAvgPlot(avg_s_rate)
%print(fullfile(figSaveDir,'blockwise_successRate'), '-dpdf', '-vector')
% success rate excluding omissions
avg_s_rate_wo_oms = cell2mat(cellfun(@(a) sum(a)/length(a), success_rate_wo_oms, 'un', 0)).*100; 
blockwiseAvgPlot(avg_s_rate_wo_oms)
ylim([20 110])
%print(fullfile(figSaveDir,'blockwise_successRate_without_omission'), '-dpdf', '-vector')
% omission rate
avg_oms_rate = cell2mat(cellfun(@(a) sum(a)/length(a), oms_rate, 'un', 0)).*100; 
blockwiseAvgPlot(avg_oms_rate)
%print(fullfile(figSaveDir,'blockwise_omissionRate'), '-dpdf', '-vector')

% repeated measures analysis of variance
rm_rez.success_rate = blockByblock_rm(avg_s_rate_wo_oms, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
success_rate_m_b_avg_C = table2cell(rm_rez.success_rate.rm_tab); 
success_rate_m_b_avg = cell2mat(success_rate_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(success_rate_m_b_avg); 
ylim([40 100])
%print(fullfile(figSaveDir,'blockwise_success_rate_m_b_combined'), '-dpdf', '-vector')

%% reach duration
med_rch_dur = cell2mat(cellfun(@nanmedian, rch_dur, 'un', 0)); 
rm_rez.rch_dur = blockByblock_rm(med_rch_dur, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
rch_dur_m_b_avg_C = table2cell(rm_rez.rch_dur.rm_tab); 
rch_dur_m_b_avg = cell2mat(rch_dur_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(rch_dur_m_b_avg); 
ylim([120 700])
print(fullfile(figSaveDir,'blockwise_rch_dur_m_b_combined'), '-dpdf', '-vector')
[avg_rch_dur, ~, sem_rch_dur] = meanstdsem(med_rch_dur); 

blockwiseAvgPlot(med_rch_dur); 
%ylim([100 550])
print(fullfile(figSaveDir,'blockwise_rch_dur_all_blocks_sessions'), '-dpdf', '-vector')

%% pull duration
med_pull_dur = cell2mat(cellfun(@nanmedian, pull_dur, 'un', 0)); 
rm_rez.pull_dur = blockByblock_rm(med_pull_dur, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
pull_dur_m_b_avg_C = table2cell(rm_rez.pull_dur.rm_tab); 
pull_dur_m_b_avg = cell2mat(pull_dur_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(pull_dur_m_b_avg); 
ylim([100 200])
print(fullfile(figSaveDir,'blockwise_pull_dur_m_b_combined'), '-dpdf', '-vector')
[avg_pull_dur, ~, sem_pull_dur] = meanstdsem(med_pull_dur); 

blockwiseAvgPlot(med_pull_dur); 
%ylim([100 550])
print(fullfile(figSaveDir,'blockwise_pull_dur_all_blocks_sessions'), '-dpdf', '-vector')

%% reach-pull duration (pulled)
med_rch_pull_dur = med_rch_dur + med_pull_dur; 
rm_rez.rch_pull_dur = blockByblock_rm(med_rch_pull_dur, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
rch_pull_dur_m_b_avg_C = table2cell(rm_rez.rch_pull_dur.rm_tab); 
rch_pull_dur_m_b_avg = cell2mat(rch_pull_dur_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(rch_pull_dur_m_b_avg); 
ylim([120 700])
print(fullfile(figSaveDir,'blockwise_reach_pull_dur_m_b_combined'), '-dpdf', '-vector')

blockwiseAvgPlot(med_rch_pull_dur); 
print(fullfile(figSaveDir,'blockwise_reach_pull_dur_all_blocks_sessions'), '-dpdf', '-vector')

%% rt
med_rt = cell2mat(cellfun(@nanmedian, rt_blk, 'un', 0)); 
rm_rez.rt = blockByblock_rm(med_rt, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
rt_m_b_avg_C = table2cell(rm_rez.rt.rm_tab); 
rt_m_b_avg = cell2mat(rt_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(rt_m_b_avg); 
print(fullfile(figSaveDir,'blockwise_rt_m_b_combined'), '-dpdf', '-vector')

blockwiseAvgPlot(med_rt); 
print(fullfile(figSaveDir,'blockwise_rt_all_blocks_sessions'), '-dpdf', '-vector')

%% reach speed
med_max_rch_spd = cell2mat(cellfun(@nanmedian, max_rch_spd, 'un', 0)); 
rm_rez.max_rch_spd = blockByblock_rm(med_max_rch_spd, filePath, 'avg_b', true, 'avg_m', true, 'b_pairs', {[1, 5], [2, 6], [3, 7], [4, 8]});
rch_spd_m_b_avg_C = table2cell(rm_rez.max_rch_spd.rm_tab); 
rch_spd_m_b_avg = cell2mat(rch_spd_m_b_avg_C(:, 2:end)); 
blockwiseAvgPlot(rch_spd_m_b_avg); 
ylim([25 40])
print(fullfile(figSaveDir,'blockwise_rch_spd_m_b_combined'), '-dpdf', '-vector')

blockwiseAvgPlot(med_max_rch_spd); 
print(fullfile(figSaveDir,'blockwise_rch_spd_all_blocks_sessions'), '-dpdf', '-vector')



%% save
save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'), 'rm_rez', '-append')
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez.mat'), 'rm_rez')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
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


function rez = blockByblock_rm(sbb0, file_path, varargin)
% sbb0: session-by-block behavioral data 

% preprocess the data matrix
if size(sbb0, 1) ~= length(file_path)
    error("Check the input data, they don't match!")
end

p = parse_input_blockByblock_rm(sbb0, file_path, varargin);

if p.Results.avg_m  % animal average
    for i = 1:length(file_path)
        wrI = strfind(file_path{i}, 'WR');
        m_name_c{i, 1} = file_path{i}(wrI:wrI+3);  
    end
    
    unique_m_name_c = unique(m_name_c); 
    sbb1 = nan(size(unique_m_name_c, 1), size(sbb0, 2)); 
    
    for j = 1:length(unique_m_name_c)
        m_I = cellfun(@(a) strcmpi(a, unique_m_name_c{j}), m_name_c); 
        sbb1(j, :) = nanmean(sbb0(m_I, :), 1);  
    end
else  % without animal average 
    sbb1 = sbb0; 
end

if p.Results.avg_b  % block average
   sbb2 = nan(size(sbb1, 1), length(p.Results.b_pairs)); 
   for j = 1:length(p.Results.b_pairs) 
       sbb2(:, j) = nanmean(sbb1(:, p.Results.b_pairs{j}), 2); 
   end
else  % without block average
    sbb2 = sbb1; 
end
    
%% run repeated measures of ANOVA
% artificial grouping to get by ranova
group_var = cell(size(sbb2, 1), 1); 
[group_var{1:ceil(size(sbb2, 1)/2)}] = deal('m1'); 
[group_var{ceil(size(sbb2, 1)/2)+1:end}] = deal('m2'); 

% variable names 
var_names = cell(1, 1+size(sbb2, 2)); 
for jj = 1:length(var_names)
    if jj == 1
        var_names{jj} = 'G'; 
    else
        var_names{jj} = strcat('B', num2str(jj-1)); 
    end
end

% build table
rm_tab = cell2table([group_var, num2cell(sbb2)], 'VariableNames', var_names); 

Block = table((1:size(sbb2,2))','VariableNames',{'Block'}); % Block
% Fit repetitive model to data
rm_design = ['B1-', strcat('B', num2str(size(sbb2,2))), ' ', '~', ' ', 'G'];  % rm design string

rm = fitrm(rm_tab, rm_design, 'WithinDesign', Block); % repeated measures model fit for RT data
ranova_rez = ranova(rm); % repeated measures ANOVA on RT (within- and within*group interaction)
multiCompBlock = multcompare(rm,'Block');

rez.rm = rm; 
rez.ranova_rez = ranova_rez; 
rez.multiCompBlock = multiCompBlock;  
rez.m_names = unique_m_name_c; 
rez.rm_tab = rm_tab; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_blockByblock_rm(sbb0, file_path, vargs)
        % sbb0: session-by-block behavioral data 
        default_avg_b = true;
        default_avg_m = true;
        default_b_pairs = {[1, 5], [2, 6], [3, 7], [4, 8]}; 
        
        p = inputParser; % create parser object
        addRequired(p, 'sbb0')
        addRequired(p, 'file_path')
        addParameter(p, 'avg_b', default_avg_b)
        addParameter(p, 'avg_m', default_avg_m)
        addParameter(p, 'b_pairs', default_b_pairs)
        
        parse(p, sbb0, file_path, vargs{:})
        
    end
end

function [rch_dur_blk, pull_dur_blk] = bByb_reach_pull_duration(ss_blk)

rch_dur_blk = nan(length(ss_blk), 1);
pull_dur_blk = nan(length(ss_blk), 1);

for jj = 1:length(ss_blk) 
    ttI = strcmpi(ss_blk(jj).trialType, 'sp'); 
    evtI = strcmpi(ss_blk(jj).evtAlign, 'rStart'); 
    timeAlign = ss_blk(jj).timeAlign;
    pullStart = ss_blk(jj).tPullStart;
    pullStop = ss_blk(jj).tPullStop;
    
    if ttI && evtI
        if ~isempty(pullStart)
            rch_dur_blk(jj, 1) = max(0, pullStart - timeAlign); 
        end
        
        if ~isempty(pullStop)
            pull_dur_blk(jj, 1) = max(0, pullStop - pullStart); 
        end
    end
end
end

