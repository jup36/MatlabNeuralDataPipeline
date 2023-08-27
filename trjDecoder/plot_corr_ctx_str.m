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
