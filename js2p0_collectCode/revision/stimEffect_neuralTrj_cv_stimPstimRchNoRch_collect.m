 filePaths = {
    '/Volumes/Extreme SSD/js2p0/WR38_052119/Matfiles', ... % Dual recording without Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR44_031020/Matfiles'};    % Dual recording with contra Cg delayed silencing (checked)

dims = 10; 
folds = 10; 

% preprocessing
%for f = 1:length(filePaths)
%    js2p0_tbytSpkHandJsPreprocess_50ms_stimPstim(filePaths{f})
%end

% train RRR with movement-potent corticostriatal activity and test them
% for stim- and pseudo-stim aligned activity
for f = 1:length(filePaths)
    filePath = GrabFiles_sort_trials('js2p0_tbytSpkHandJsTrjBin_50ms_stimPstimRchNoRch_', 0, filePaths(f));
    [~, rrrRezR2] = stimEffect_neuralTrj_cv_stim_pstim_rch_noRch(filePath{1}, 20, 10); % args: dimensions and folds

    % organize r2 results
    mR2.noStim_tbyt{f} = rrrRezR2.mR2_tbyt;
    mR2.stimRch_tbyt{f}   = rrrRezR2.mR2_tbytStimRch;
    mR2.stimNoRch_tbyt{f}   = rrrRezR2.mR2_tbytStimNoRch;
    mR2.pStimRch_tbyt{f}   = rrrRezR2.mR2_tbytPstimRch;
    mR2.pStimNoRch_tbyt{f}   = rrrRezR2.mR2_tbytPstimNoRch;   
  
    mR2.noStim{f} = rrrRezR2.mR2;
    mR2.stimRch{f}   = rrrRezR2.mR2_stimRch;
    mR2.stimNoRch{f}   = rrrRezR2.mR2_stimNoRch;
    mR2.pStimRch{f}   = rrrRezR2.mR2_pStimRch;
    mR2.pStimNoRch{f}   = rrrRezR2.mR2_pStimNoRch;   

    medR2.noStim_tbyt{f} = rrrRezR2.medR2_tbyt;
    medR2.stimRch_tbyt{f}   = rrrRezR2.medR2_tbytStimRch;
    medR2.stimNoRch_tbyt{f}   = rrrRezR2.medR2_tbytStimNoRch;
    medR2.pStimRch_tbyt{f}   = rrrRezR2.medR2_tbytPstimRch;
    medR2.pStimNoRch_tbyt{f}   = rrrRezR2.medR2_tbytPstimNoRch;   
  
    medR2.noStim{f} = rrrRezR2.medR2;
    medR2.stimRch{f}   = rrrRezR2.medR2_stimRch;
    medR2.stimNoRch{f}   = rrrRezR2.medR2_stimNoRch;
    medR2.pStimRch{f}   = rrrRezR2.medR2_pStimRch;
    medR2.pStimNoRch{f}   = rrrRezR2.medR2_pStimNoRch;   

    fprintf("Completed file #%d\n", f)
end

save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', ['stimEffect_neuralTrj_cv_collect_stimPstimRchNoRch', sprintf('_Dims%d', dims), sprintf('_Folds%d', folds)]), 'mR2*', 'medR2*')

%% Stats
% prepare session-by-epoch data
mR2_noStim_pStim_Stim_prep = [cell2mat(mR2.noStim)', cell2mat(mR2.pStim)', cell2mat(mR2.stim)', cell2mat(mR2.prep)']; 
mR2_tbyt_noStim_pStim_Stim_prep = [cell2mat(mR2.noStim_tbyt)', cell2mat(mR2.pStim_tbyt)', cell2mat(mR2.stim_tbyt)', cell2mat(mR2.prep_tbyt)']; 
medR2_noStim_pStim_Stim_prep = [cell2mat(medR2.noStim)', cell2mat(medR2.pStim)', cell2mat(medR2.stim)', cell2mat(medR2.prep)']; 
medR2_tbyt_noStim_pStim_Stim_prep = [cell2mat(medR2.noStim_tbyt)', cell2mat(medR2.pStim_tbyt)', cell2mat(medR2.stim_tbyt)', cell2mat(medR2.prep_tbyt)']; 

% run stats
stat.mR2 = onewayANOVA_multiCompare(mR2_noStim_pStim_Stim_prep); 
stat.mR2_tbyt = onewayANOVA_multiCompare(mR2_tbyt_noStim_pStim_Stim_prep); 
stat.medR2 = onewayANOVA_multiCompare(medR2_noStim_pStim_Stim_prep); 
stat.medR2_tbyt = onewayANOVA_multiCompare(medR2_tbyt_noStim_pStim_Stim_prep); 

save(fullfile('/Volumes/Extreme SSD/js2p0/collectData', ['stimEffect_neuralTrj_cv_collect_stimPstimPrepExt', sprintf('_Dims%d', dims), sprintf('_Folds%d', folds)]), 'stat', '-append')

%% plot
fig = scatter_row_by_row([cell2mat(mR2_pStim)', cell2mat(mR2)', cell2mat(mR2_stim)'], [0, 0.4]); 
print(fig, fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', ...
    ['stimEffect_neuralTrj_cv_collect_rez_mR2', sprintf('_Dims%d', dims), sprintf('_Folds%d', folds)]), ...
    '-vector', '-dpdf');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig = scatter_row_by_row(collect_mat, ylimits)

fig = figure; hold on; 
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


function rez = blockByblock_rm_simple(sbb0)
% sbb0: session-by-block data 
    
%% run repeated measures of ANOVA
% artificial grouping to get by ranova
group_var = cell(size(sbb0, 1), 1); 
[group_var{1:ceil(size(sbb0, 1)/2)}] = deal('m1'); 
[group_var{ceil(size(sbb0, 1)/2)+1:end}] = deal('m2'); 

% variable names 
var_names = cell(1, 1+size(sbb0, 2)); 
for jj = 1:length(var_names)
    if jj == 1
        var_names{jj} = 'G'; 
    else
        var_names{jj} = strcat('B', num2str(jj-1)); 
    end
end

% build table
rm_tab = cell2table([group_var, num2cell(sbb0)], 'VariableNames', var_names); 

Block = table((1:size(sbb0,2))','VariableNames',{'Block'}); % Block
% Fit repetitive model to data
rm_design = ['B1-', strcat('B', num2str(size(sbb0,2))), ' ', '~', ' ', 'G'];  % rm design string

rm = fitrm(rm_tab, rm_design, 'WithinDesign', Block); % repeated measures model fit for RT data
ranova_rez = ranova(rm); % repeated measures ANOVA on RT (within- and within*group interaction)
multiCompBlock = multcompare(rm,'Block');

rez.rm = rm; 
rez.ranova_rez = ranova_rez; 
rez.multiCompBlock = multiCompBlock;  
rez.rm_tab = rm_tab; 

end

function one = onewayANOVA_multiCompare(dat)
%This function runs one-way anova and posthoc pairwise test with correction
% for multiple comparisons.
% Input: dat is n-by-p matrix (n: data points, p: groups).
% Output: one (output structure). Note that the multiple comparison result
% is contained in one.c, see below for interpretation of c.
%For example, suppose one row contains the following entries.
% 2.0000  5.0000  1.9442  8.2206  14.4971 0.0432
% These numbers indicate that the mean of group 2 minus the mean of group 5 is estimated to be 8.2206,
% and a 95% confidence interval for the true difference of the means is [1.9442, 14.4971].
% The p-value for the corresponding hypothesis test that the difference of the means of groups 2 and 5
% is significantly different from zero is 0.0432.

% One-way ANOVA
[one.p, one.tbl, one.stats] = anova1(dat);
[one.c, ~, ~, ~] = multcompare(one.stats, "CriticalValueType", "tukey-kramer");
end

function fig = meanSemErrorbar(mean_values, sem_values)

fig = figure; 
% Assuming mean_values and sem_values are already defined as 1x20 vectors

% Generate a vector for the x-axis
x_values = 1:length(mean_values); 

% Plot filled circles for mean values
scatter(x_values, mean_values, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');

% Hold on to the current figure
hold on;

% Add error bars for standard error of the mean
errorbar(x_values, mean_values, sem_values, 'LineStyle', 'none', 'Color', 'r', 'CapSize', 10);

% Label the axes
xlabel('Dims');
ylabel('Mean Value');

% Add a title and a legend
%title('Mean Values with SEM');
%legend('Mean values', 'SEM', 'Location', 'best');

% Hold off to finish the plotting
hold off;

end







