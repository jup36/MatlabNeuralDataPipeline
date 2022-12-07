function plot_corr_ctx_str(corr_cell, ylimits)

ctxI = cell2mat(cellfun(@(a) strcmp(a, 'ctx'), corr_cell(1,:), 'un', 0)); 
strI = cell2mat(cellfun(@(a) strcmp(a, 'str'), corr_cell(1,:), 'un', 0)); 
ctx_strI = cell2mat(cellfun(@(a) strcmp(a, 'ctx_str'), corr_cell(1,:), 'un', 0)); 

ctx_corr = corr_cell(2:end, ctxI); 
str_corr = corr_cell(2:end, strI); 
ctx_str_corr = corr_cell(2:end, ctx_strI); 

collect_corr_C = [ctx_corr, str_corr, ctx_str_corr]; 
empties = cell2mat(cellfun(@isempty, collect_corr_C, 'un', 0));

[collect_corr_C{empties}] = deal(NaN); 
collect_corr = cell2mat(cellfun(@(a) nanmean(a), collect_corr_C, 'un', 0)); 

scatter_row_by_row(collect_corr, ylimits)

end
