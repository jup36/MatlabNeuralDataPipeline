function plot_corr_ctx_str_cg(corr_cell, ylimits)

ctxI = cell2mat(cellfun(@(a) strcmp(a, 'ctx1'), corr_cell(1,:), 'un', 0)); 
strI = cell2mat(cellfun(@(a) strcmp(a, 'str1'), corr_cell(1,:), 'un', 0)); 
cgI = cell2mat(cellfun(@(a) strcmp(a, 'cg'), corr_cell(1,:), 'un', 0)); 

ctx_corr = corr_cell(2:end, ctxI); 
ctx_corr_I = ~cell2mat(cellfun(@isempty, ctx_corr, 'un', 0)); 
str_corr = corr_cell(2:end, strI); 
str_corr_I = ~cell2mat(cellfun(@isempty, str_corr, 'un', 0)); 
cg_corr = corr_cell(2:end, cgI); 
cg_corr_I = ~cell2mat(cellfun(@isempty, cg_corr, 'un', 0)); 

collect_corr_C = [ctx_corr, str_corr, cg_corr]; 
empties = cell2mat(cellfun(@isempty, collect_corr_C, 'un', 0));

[collect_corr_C{empties}] = deal(NaN); 
collect_corr = cell2mat(cellfun(@(a) nanmean(a), collect_corr_C, 'un', 0)); 

scatter_row_by_row(collect_corr, ylimits)

end