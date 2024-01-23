filePath = {'/Volumes/Extreme SSD/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22) % GOOD TO BE USED
            '/Volumes/Extreme SSD/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked            % GOOD TO BE USED
            '/Volumes/Extreme SSD/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked          % BAD (Lots of Matrix singular or bad-scaled warnings)
            '/Volumes/Extreme SSD/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked            % There was an issue with 'corr' after running iterations that needs to be revisited!
            '/Volumes/Extreme SSD/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked            % GOOD
            '/Volumes/Extreme SSD/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked       
            '/Volumes/Extreme SSD/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
            '/Volumes/Extreme SSD/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/Extreme SSD/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked      
            '/Volumes/Extreme SSD/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/Extreme SSD/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
            '/Volumes/Extreme SSD/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked % GOOD
        
pos_corr = cell(length(filePath)+1, 6+1);
pos_corr(1,:) = deal({'cat', 'ctx', 'ctx_mismatch', 'str', 'str_mismatch', 'cg', 'cg_mismatch'});

pos_r2 = cell(length(filePath)+1, 6+1);
pos_r2(1,:) = deal({'cat', 'ctx', 'ctx_mismatch', 'str', 'str_mismatch', 'cg', 'cg_mismatch'});

for f = 1:length(filePath)
    fileInfo = dir(fullfile(filePath{f}, 'Matfiles','rezKFdecodeHTrjCtxStrPosVel_weight_mismatch_*'));
    if ~isempty(fileInfo)
        %% load reach data
        load(fullfile(fileInfo.folder, fileInfo.name), 'corrRez', 'r2Rez');
        % reach position
        name_loc = strfind(fileInfo.name, 'WR');
        save_name = fileInfo.name(name_loc:name_loc+10);
        pos_corr = organize_corr_weight_switch(corrRez, pos_corr, f+1, save_name);
        pos_r2 = organize_r2_weight_switch(r2Rez, pos_r2, f+1, save_name);
    end
    fprintf('processed file #%d\n', f)
end

%% ctx weight match vs. mismatch
r2.ctx_match_all = cell2mat(pos_r2(2:end, 2)); 
r2.ctx_match = r2.ctx_match_all(:, end); 
r2.ctx_match(r2.ctx_match<0) = eps; 

r2.ctx_mismatch_all = cell2mat(pos_r2(2:end, 3)); 
r2.ctx_mismatch = r2.ctx_mismatch_all(:, end); 
r2.ctx_mismatch(r2.ctx_mismatch<0) = eps; 

[~, r2.ctx_p, ~, r2.ctx_stats] = ttest(r2.ctx_match, r2.ctx_mismatch); 
scatter_row_by_row_2col([r2.ctx_match, r2.ctx_mismatch], [0 0.6])

%% str weight match vs. mismatch
r2.str_match_all = cell2mat(pos_r2(2:end, 4)); 
r2.str_match = r2.str_match_all(:, end); 
r2.str_match(r2.str_match<0) = eps; 

r2.str_mismatch_all = cell2mat(pos_r2(2:end, 5)); 
r2.str_mismatch = r2.str_mismatch_all(:, end); 
r2.str_mismatch(r2.str_mismatch<0) = eps; 

[~, r2.str_p, ~, r2.str_stats] = ttest(r2.str_match, r2.str_mismatch); 
scatter_row_by_row_2col([r2.str_match, r2.str_mismatch], [0 0.6])

%% cg weight match vs. mismatch
r2.cg_match_all = cell2mat(pos_r2(2:end, 6)); 
r2.cg_match = r2.cg_match_all(:, end); 
r2.cg_match(r2.cg_match<0) = eps; 

r2.cg_mismatch_all = cell2mat(pos_r2([3, 5, 6, 13, 14], 7)); 
r2.cg_mismatch = r2.cg_mismatch_all(:, end); 
r2.cg_mismatch(r2.cg_mismatch<0) = eps; 

[~, r2.cg_p, ~, r2.cg_stats] = ttest(r2.cg_match, r2.cg_mismatch); 
scatter_row_by_row_2col([r2.cg_match, r2.cg_mismatch], [0 0.6])

%% ctx, str, cg weight match vs. mismatch
[~, r2.ctx_str_cg_p, ~, r2.ctx_str_cg_stats] = ttest([r2.ctx_match; r2.str_match; r2.cg_match], [r2.ctx_mismatch; r2.str_mismatch; r2.cg_mismatch]); 
scatter_row_by_row_2col([[r2.ctx_match; r2.str_match; r2.cg_match], [r2.ctx_mismatch; r2.str_mismatch; r2.cg_mismatch]], [0 0.6])
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_cg_weight_match_mismatch'), '-dpdf', '-painters')

mean([r2.ctx_mismatch; r2.str_mismatch; r2.cg_mismatch]-[r2.ctx_match; r2.str_match; r2.cg_match])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrRez_collect = organize_corr_weight_switch(corrRez, corrRez_collect, row_I, saveName)

fields = fieldnames(corrRez); 

corrRez_collect{row_I, 1} = saveName;
for j = 1:length(fields)
    field_name = fields{j}; 
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), corrRez_collect(1,:), 'un', 0)); 
    if isstruct(corrRez.(field_name))
        corr_xyz = diag(corrRez.(field_name).all_sm)'; 
    else
        corr_xyz = []; 
    end
    corrRez_collect{row_I, col_I} = corr_xyz; 
end

end

function r2Rez_collect = organize_r2_weight_switch(r2Rez, r2Rez_collect, row_I, saveName)

fields = fieldnames(r2Rez); 

r2Rez_collect{row_I, 1} = saveName;
for j = 1:length(fields)
    field_name = fields{j}; 
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), r2Rez_collect(1,:), 'un', 0)); 
    if isstruct(r2Rez.(field_name))
        r2_xyz_all = [r2Rez.(field_name).all, r2Rez.(field_name).overall]; 
    else
        r2_xyz_all = []; 
    end
    r2Rez_collect{row_I, col_I} = r2_xyz_all; 
end

end

function plot_within_reach_vs_pull(reach_cell, pull_cell, region, ylimits)
regionI = cell2mat(cellfun(@(a) strcmp(a, region), reach_cell(1,:), 'un', 0));

reach_region = reach_cell(2:end, regionI); 
pull_region = pull_cell(2:end, regionI); 

collect_C = [reach_region, pull_region]; 
empties = cell2mat(cellfun(@isempty, collect_C, 'un', 0));

[collect_C{empties}] = deal(NaN); 
collect_region = cell2mat(cellfun(@(a) a(end), collect_C, 'un', 0)); 

scatter_row_by_row(collect_region, ylimits)

end


function scatter_row_by_row_2col(collect_mat, ylimits)

figure; hold on; 
for r = 1:size(collect_mat, 1) 
    x = rand(1, size(collect_mat, 2)).*0.1-0.05 + [1:size(collect_mat, 2)]; 
    plot(x, collect_mat(r, :), 'k:')
    scatter(x, collect_mat(r, :), 50)
end
clearvars r

ylim(ylimits)

xlim([.5 2.5])
set(gca, 'TickDir', 'out')

avg_collect_mat = nanmean(collect_mat(sum(collect_mat>0, 2)==size(collect_mat, 2), :), 1); 
for rr = 1:size(collect_mat, 2)
    plot([rr-.2, rr+.2], avg_collect_mat(rr).*ones(1, 2), 'k', 'LineWidth', 2.5)
end

end

function collect_R2 = plot_r2_ctx_str_cg(r2_cell, ylimits)

ctxI = cell2mat(cellfun(@(a) strcmp(a, 'ctx1'), r2_cell(1,:), 'un', 0)); 
strI = cell2mat(cellfun(@(a) strcmp(a, 'str1'), r2_cell(1,:), 'un', 0)); 
cgI = cell2mat(cellfun(@(a) strcmp(a, 'cg'), r2_cell(1,:), 'un', 0)); 

ctx_R2 = r2_cell(2:end, ctxI); 
ctx_R2_I = ~cell2mat(cellfun(@isempty, ctx_R2, 'un', 0)); 
str_R2 = r2_cell(2:end, strI); 
str_R2_I = ~cell2mat(cellfun(@isempty, str_R2, 'un', 0)); 
cg_R2 = r2_cell(2:end, cgI); 
cg_R2_I = ~cell2mat(cellfun(@isempty, cg_R2, 'un', 0)); 

collect_R2_C = [ctx_R2, str_R2, cg_R2]; 
empties = cell2mat(cellfun(@isempty, collect_R2_C, 'un', 0));

[collect_R2_C{empties}] = deal(NaN); 
collect_R2 = cell2mat(cellfun(@(a) a(end), collect_R2_C, 'un', 0)); 

scatter_row_by_row(collect_R2, ylimits)

end







