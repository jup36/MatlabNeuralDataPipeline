filePath = {'/Volumes/Extreme SSD/js2p0/WR37_022119', ... % Dual recording without silencing, Trj checked (07/04/22) % GOOD TO BE USED
    '/Volumes/Extreme SSD/js2p0/WR37_022219', ... % Dual recording without silencing, Trj checked            % BAD (Lots of Matrix singular or bad-scaled warnings)
    '/Volumes/Extreme SSD/js2p0/WR37_022619', ... % Cg recording contra-Cg silencing, Trj checked            % GOOD TO BE USED
    '/Volumes/Extreme SSD/js2p0/WR37_022719', ... % B Only contra-Cg delayed silencing, Trj checked          % BAD (Lots of Matrix singular or bad-scaled warnings)
    '/Volumes/Extreme SSD/js2p0/WR38_052119', ... % Dual recording without silencing, Trj checked            % There was an issue with 'corr' after running iterations that needs to be revisited!
    '/Volumes/Extreme SSD/js2p0/WR38_052219', ... % Dual recording without silencing, Trj checked            % GOOD
    '/Volumes/Extreme SSD/js2p0/WR38_052319', ... % Cg recording contra-Cg silencing, Trj checked            % BAD (Crashed)
    '/Volumes/Extreme SSD/js2p0/WR38_052419', ... % Corticostriatal recording M1 silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_091019', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_091119', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR39_100219', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
    '/Volumes/Extreme SSD/js2p0/WR39_100319', ... % B Only contra-Cg delayed silencing, Trj checked
    '/Volumes/Extreme SSD/js2p0/WR40_081919', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
    '/Volumes/Extreme SSD/js2p0/WR40_082019', ... % Dual recording with contra Cg silencing, Trj checked     % GOOD
    '/Volumes/Extreme SSD/js2p0/WR44_031020'};    % Dual recording with contra Cg delayed silencing, Trj checked % GOOD

reach_pos_corr = cell(length(filePath)+1, 7+1);
reach_pos_corr(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

reach_pos_r2 = cell(length(filePath)+1, 7+1);
reach_pos_r2(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

reach_vel_corr = cell(length(filePath)+1, 7+1);
reach_vel_corr(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

reach_vel_r2 = cell(length(filePath)+1, 7+1);
reach_vel_r2(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

pull_pos_corr = cell(length(filePath)+1, 7+1);
pull_pos_corr(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

pull_pos_r2 = cell(length(filePath)+1, 7+1);
pull_pos_r2(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

pull_vel_corr = cell(length(filePath)+1, 7+1);
pull_vel_corr(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

pull_vel_r2 = cell(length(filePath)+1, 7+1);
pull_vel_r2(1,:) = deal({'cat', 'ctx', 'ctx1', 'str', 'str1', 'cg', 'ctx_str', 'ctx_str_cg'});

for f = 1:length(filePath)
    fileInfo = dir(fullfile(filePath{f}, 'Matfiles','rezKFdecodeHTrjCtxStrPosVel_reach_new_*'));
    if ~isempty(fileInfo)
        %% load reach data
        load(fullfile(fileInfo.folder, fileInfo.name), 'corrRez_reach', 'r2Rez_reach');
        % reach position
        name_loc = strfind(fileInfo.name, 'WR');
        save_name = fileInfo.name(name_loc:name_loc+10);
        reach_pos_corr = organize_corr(corrRez_reach, reach_pos_corr, f+1, 7, save_name);
        reach_pos_r2 = organize_r2(r2Rez_reach, reach_pos_r2, f+1, 7, save_name);
        % reach velocity
        reach_vel_corr = organize_corr(corrRez_reach, reach_vel_corr, f+1, 8, save_name);
        reach_vel_r2 = organize_r2(r2Rez_reach, reach_vel_r2, f+1, 8, save_name);
    end
    
    %% load pull data
    fileInfo = dir(fullfile(filePath{f}, 'Matfiles','rezKFdecodeHTrjCtxStrPosVel_pull_new_*'));
    if ~isempty(fileInfo)
        load(fullfile(fileInfo.folder, fileInfo.name), 'corrRez_pull', 'r2Rez_pull');
        name_loc = strfind(fileInfo.name, 'WR');
        save_name = fileInfo.name(name_loc:name_loc+10);
        % pull position
        pull_pos_corr = organize_corr(corrRez_pull, pull_pos_corr, f+1, 7, save_name);
        pull_pos_r2 = organize_r2(r2Rez_pull, pull_pos_r2, f+1, 7, save_name);
        % pull velocity
        pull_vel_corr = organize_corr(corrRez_pull, pull_pos_corr, f+1, 8, save_name);
        pull_vel_r2 = organize_r2(r2Rez_pull, pull_pos_r2, f+1, 8, save_name);
    end
    fprintf('processed file #%d\n', f)
end

%% reach position overall R2, CTX-STR
r2_ctx_str_reach = plot_r2_ctx_str(reach_pos_r2, [.2 .8], 4);  % plot ctx-str R2 reach
[anova_rez.r2_ctx_str_reach_dat, anova_rez.r2_ctx_str_reach_label] = organize_anova_cell(r2_ctx_str_reach); 
[anova_rez.r2_ctx_str_reach.p, anova_rez.r2_ctx_str_reach.tbl, anova_rez.r2_ctx_str_reach.stats] = anova1(anova_rez.r2_ctx_str_reach_dat, anova_rez.r2_ctx_str_reach_label); 
%print(fullfile('/Users/parkj/Desktop/js2p0_temp_ncm', 'reach_pos_r2_ctx_str'), '-dpdf', '-vector')

%% reach X position R2, CTX-STR
r2_ctx_str_reach_X = plot_r2_ctx_str_x(reach_pos_r2, [.2 .8]);  % plot ctx-str R2 reach
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_X'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_reach_X.p, anova_rez.r2_ctx_str_reach_X.tbl, anova_rez.r2_ctx_str_reach_X.stats] = anova1(r2_ctx_str_reach_X); 

%% reach Y position R2, CTY-STR
r2_ctx_str_reach_Y = plot_r2_ctx_str_y(reach_pos_r2, [.2 .8]);  % plot ctx-str R2 reach
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_Y'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_reach_Y.p, anova_rez.r2_ctx_str_reach_Y.tbl, anova_rez.r2_ctx_str_reach_Y.stats] = anova1(r2_ctx_str_reach_Y); 

%% reach Z position R2, CTZ-STR
r2_ctx_str_reach_Z = plot_r2_ctx_str_z(reach_pos_r2, [.2 .8]);  % plot ctx-str R2 reach
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_Z'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_reach_Z.p, anova_rez.r2_ctx_str_reach_Z.tbl, anova_rez.r2_ctx_str_reach_Z.stats] = anova1(r2_ctx_str_reach_Z); 

%% pull overall R2, CTX-STR
r2_ctx_str_pull = plot_r2_ctx_str(pull_pos_r2, [.2 .8], 4);  % plot ctx-str R2 pull
[anova_rez.r2_ctx_str_pull.p, anova_rez.r2_ctx_str_pull.tbl, anova_rez.r2_ctx_str_pull.stats] = anova1(r2_ctx_str_pull); 

%plot_corr_ctx_str(reach_pos_corr, [.2 .8]) % plot ctx-str R2 reach
%print(fullfile('/Users/parkj/Desktop/js2p0_temp_ncm', 'reach_pos_corr_ctx_str'), '-dpdf', '-vector')

%% pull X position R2, CTX-STR
r2_ctx_str_pull_X = plot_r2_ctx_str_x(pull_pos_r2, [.2 .8]);  % plot ctx-str R2 pull
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_X'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_pull_X.p, anova_rez.r2_ctx_str_pull_X.tbl, anova_rez.r2_ctx_str_pull_X.stats] = anova1(r2_ctx_str_pull_X); 

%% pull Y position R2, CTY-STR
r2_ctx_str_pull_Y = plot_r2_ctx_str_y(pull_pos_r2, [.2 .8]);  % plot ctx-str R2 pull
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_Y'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_pull_Y.p, anova_rez.r2_ctx_str_pull_Y.tbl, anova_rez.r2_ctx_str_pull_Y.stats] = anova1(r2_ctx_str_pull_Y); 

%% pull Z position R2, CTZ-STR
r2_ctx_str_pull_Z = plot_r2_ctx_str_z(pull_pos_r2, [.2 .8]);  % plot ctx-str R2 pull
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_Z'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_pull_Z.p, anova_rez.r2_ctx_str_pull_Z.tbl, anova_rez.r2_ctx_str_pull_Z.stats] = anova1(r2_ctx_str_pull_Z); 

%% reach overall R2, CTX-STR-Cg
r2_ctx_str_cg_reach = plot_r2_ctx_str_cg(reach_pos_r2, [.1 .8], 4);  % plot ctx-str-cg R2 reach overall
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_cg_XYZ'), '-dpdf', '-vector')

[anova_rez.r2_ctx_str_cg_reach.p, anova_rez.r2_ctx_str_cg_reach.tbl, anova_rez.r2_ctx_str_cg_reach.stats] = anova1(r2_ctx_str_cg_reach); 
anova_rez.r2_ctx_str_cg_reach.multicompare = multcompare(anova_rez.r2_ctx_str_cg_reach.stats, "ctype", "bonferroni");
anova_rez.r2_ctx_str_cg_reach.tbl = array2table(anova_rez.r2_ctx_str_cg_reach.multicompare); 

%% reach X position R2, CTX-STR-Cg
r2_ctx_str_cg_reach_X = plot_r2_ctx_str_cg_X(reach_pos_r2, [0 .8]);  % plot ctx-str-cg R2 reach X position
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_cg_X'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_cg_reach_X.p, anova_rez.r2_ctx_str_cg_reach_X.tbl, anova_rez.r2_ctx_str_cg_reach_X.stats] = anova1(r2_ctx_str_cg_reach_X); 

%% reach Y position R2, CTY-STR-Cg
r2_ctx_str_cg_reach_Y = plot_r2_ctx_str_cg_Y(reach_pos_r2, [0 .8]);  % plot ctx-str-cg R2 reach Y position
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_cg_Y'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_cg_reach_Y.p, anova_rez.r2_ctx_str_cg_reach_Y.tbl, anova_rez.r2_ctx_str_cg_reach_Y.stats] = anova1(r2_ctx_str_cg_reach_Y);  

%% reach Z position R2, CTZ-STR-Cg
r2_ctx_str_cg_reach_Z = plot_r2_ctx_str_cg_Z(reach_pos_r2, [0 .8]);  % plot ctx-str-cg R2 reach Z position
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'reach_pos_r2_ctx_str_cg_Z'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_cg_reach_Z.p, anova_rez.r2_ctx_str_cg_reach_Z.tbl, anova_rez.r2_ctx_str_cg_reach_Z.stats] = anova1(r2_ctx_str_cg_reach_Z);  

%% pull X position R2, CTX-STR-Cg
r2_ctx_str_cg_pull_X = plot_r2_ctx_str_cg_pull_X(pull_pos_r2, [0 .8]);  % plot ctx-str-cg R2 pull X position
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_cg_X'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_cg_pull_X.p, anova_rez.r2_ctx_str_cg_pull_X.tbl, anova_rez.r2_ctx_str_cg_pull_X.stats] = anova1(r2_ctx_str_cg_pull_X); 

%% pull Y position R2, CTY-STR-Cg
r2_ctx_str_cg_pull_Y = plot_r2_ctx_str_cg_pull_Y(pull_pos_r2, [0 .9]);  % plot ctx-str-cg R2 pull Y position
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_cg_Y'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_cg_pull_Y.p, anova_rez.r2_ctx_str_cg_pull_Y.tbl, anova_rez.r2_ctx_str_cg_pull_Y.stats] = anova1(r2_ctx_str_cg_pull_Y);  

%% pull Z position R2, CTZ-STR-Cg
r2_ctx_str_cg_pull_Z = plot_r2_ctx_str_cg_pull_Z(pull_pos_r2, [0 .8]);  % plot ctx-str-cg R2 pull Z position
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/KF_decoder', 'pull_pos_r2_ctx_str_cg_Z'), '-dpdf', '-vector')
[anova_rez.r2_ctx_str_cg_pull_Z.p, anova_rez.r2_ctx_str_cg_pull_Z.tbl, anova_rez.r2_ctx_str_cg_pull_Z.stats] = anova1(r2_ctx_str_cg_pull_Z);  

%% pull overall R2, CTX-STR-Cg
r2_ctx_str_cg_pull = plot_r2_ctx_str_cg(pull_pos_r2, [0 .8]);  % plot ctx-str-cg R2 pull
%print(fullfile('/Users/parkj/Desktop/js2p0_temp_ncm', 'pull_pos_r2_ctx_str_cg'), '-dpdf', '-vector')

nan_I = (isnan(r2_ctx_str_cg_pull)) | (isnan(r2_ctx_str_cg_reach)); 
r2_ctx_str_cg_reach(nan_I) = NaN;  
r2_ctx_str_cg_pull(nan_I) = NaN; 

[anova_rez.r2_ctx_str_cg_pull.p, anova_rez.r2_ctx_str_cg_pull.tbl, anova_rez.r2_ctx_str_cg_pull.stats] = anova1(r2_ctx_str_cg_pull); 
anova_rez.r2_ctx_str_cg_pull.multicompare = multcompare(anova_rez.r2_ctx_str_cg_pull.stats, "ctype", "bonferroni");
anova_rez.r2_ctx_str_cg_pull.tbl = array2table(anova_rez.r2_ctx_str_cg_pull.multicompare); 

%plot_corr_ctx_str_cg(reach_pos_corr, [0.4 0.9])
%print(fullfile('/Users/parkj/Desktop/js2p0_temp_ncm', 'reach_pos_corr_ctx_str_cg'), '-dpdf', '-vector')

%% Comparing reach vs pull decoding accuracies within each region
plot_within_reach_vs_pull(reach_pos_r2, pull_pos_r2, 'ctx', [.2 .8])

plot_within_reach_vs_pull(reach_pos_r2, pull_pos_r2, 'str', [.2 .8])

plot_within_reach_vs_pull(reach_pos_r2, pull_pos_r2, 'cg', [.1 .6])

save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'kfDecode_r2_rez'), 'anova_rez', 'reach_pos_r2', 'pull_pos_r2')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrRez_collect = organize_corr(corrRez, corrRez_collect, row_I, cell_I, saveName)

fields = fieldnames(corrRez); 

corrRez_collect{row_I, 1} = saveName;
for j = 1:length(fields)
    field_name = fields{j}; 
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), corrRez_collect(1,:), 'un', 0)); 
    if isstruct(corrRez.(field_name){1,cell_I})
        corr_xyz = diag(corrRez.(field_name){1,cell_I}.all_sm)'; 
    else
        corr_xyz = []; 
    end
    corrRez_collect{row_I, col_I} = corr_xyz; 
end

end

function r2Rez_collect = organize_r2(r2Rez, r2Rez_collect, row_I, cell_I, saveName)

fields = fieldnames(r2Rez);

r2Rez_collect{row_I, 1} = saveName;
for j = 1:length(fields)
    field_name = fields{j};
    
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), r2Rez_collect(1,:), 'un', 0));
    if isstruct(r2Rez.(field_name){1,cell_I})
        r2_xyz_all = [r2Rez.(field_name){1,cell_I}.all, r2Rez.(field_name){1,cell_I}.overall];
    else
        r2_xyz_all = [];
    end
    r2Rez_collect{row_I, col_I} = r2_xyz_all;
end

end

function scatter_row_by_row(collect_mat, ylimits)

figure; hold on; 
for r = 1:size(collect_mat, 1) 
    x = rand(1, size(collect_mat, 2)).*0.1-0.05 + (1:size(collect_mat, 2)); 
    plot(x, collect_mat(r, :), 'k:')
    scatter(x, collect_mat(r, :), 50)
end
clearvars r


ylim(ylimits)

xlim([.5 3.5])
set(gca, 'TickDir', 'out')

avg_collect_mat = nanmean(collect_mat(sum(collect_mat>0, 2)==size(collect_mat, 2), :), 1); 
for rr = 1:size(collect_mat, 2)
    plot([rr-.2, rr+.2], avg_collect_mat(rr).*ones(1, 2), 'k', 'LineWidth', 2.5)
end

end

function scatter_row_by_row_cell(collect_cell, ylimits)

figure; hold on; 
for r = 1:size(collect_cell, 1) 
    x = rand(1, size(collect_cell, 2)).*0.1-0.05 + (1:size(collect_cell, 2)); 
    plot(x, cell2mat(collect_cell(r, :)), 'k:')
    scatter(x, cell2mat(collect_cell(r, :)), 50)
end
clearvars r


ylim(ylimits)

xlim([.5 3.5])
set(gca, 'TickDir', 'out')

avg_collect_cell = nanmean(cell2mat(collect_cell)); 
for rr = 1:size(collect_cell, 2)
    plot([rr-.2, rr+.2], avg_collect_cell(rr).*ones(1, 2), 'k', 'LineWidth', 2.5)
end

end

function [anova_dat, anova_label] = organize_anova_cell(outcell)

anova_dat = []; 
anova_label = {}; 
for jj = 1:size(outcell, 2)
    % get data
    dat_col = cell2mat(outcell(:, jj)); 
    dat_col = dat_col(~isnan(dat_col));
    anova_dat = [anova_dat; dat_col];
    
    % get label
    label_col = cell(length(dat_col), 1); 
    label_col(:) = {num2str(jj)}; 
    anova_label = [anova_label; label_col]; 
end

end

function outcell = plot_corr_ctx_str(corr_cell, ylimits, col)

ctxI = cell2mat(cellfun(@(a) strcmp(a, 'ctx'), corr_cell(1,:), 'un', 0)); 
strI = cell2mat(cellfun(@(a) strcmp(a, 'str'), corr_cell(1,:), 'un', 0)); 
ctx_strI = cell2mat(cellfun(@(a) strcmp(a, 'ctx_str'), corr_cell(1,:), 'un', 0)); 

ctx_corr = corr_cell(2:end, ctxI); 
str_corr = corr_cell(2:end, strI); 
ctx_str_corr = corr_cell(2:end, ctx_strI); 

collect_corr_C = [ctx_corr, str_corr, ctx_str_corr]; 
outcell = cell(size(collect_corr_C, 1), size(collect_corr_C, 2)); 
non_nans = cell2mat(cellfun(@(a) ~isempty(sum(a, 2)), collect_corr_C, 'un', 0)); 
[outcell{~non_nans}] = deal(NaN); 
collect_corr_val = cellfun(@(a) a(col), collect_corr_C(non_nans), 'un', 0); 

outcell(non_nans) = collect_corr_val; 
outcell(cell2mat(cellfun(@(a) a<0, outcell, 'un', 0))) = {NaN}; % take negative values as NaN

scatter_row_by_row_cell(outcell, ylimits)

end








