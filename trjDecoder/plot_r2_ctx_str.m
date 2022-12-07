function collect_R2_val = plot_r2_ctx_str(r2_cell, ylimits)

ctxI = cell2mat(cellfun(@(a) strcmp(a, 'ctx'), r2_cell(1,:), 'un', 0)); 
strI = cell2mat(cellfun(@(a) strcmp(a, 'str'), r2_cell(1,:), 'un', 0)); 
ctx_strI = cell2mat(cellfun(@(a) strcmp(a, 'ctx_str'), r2_cell(1,:), 'un', 0)); 

ctx_R2 = r2_cell(2:end, ctxI); 
str_R2 = r2_cell(2:end, strI); 
ctx_str_R2 = r2_cell(2:end, ctx_strI); 

collect_R2_C = [ctx_R2, str_R2, ctx_str_R2]; 
empties = cell2mat(cellfun(@isempty, collect_R2_C, 'un', 0));

[collect_R2_C{empties}] = deal(NaN); 
collect_R2 = cell2mat(cellfun(@(a) a(end), collect_R2_C, 'un', 0)); 
non_nans = ~isnan(sum(collect_R2, 2)); 
non_negatives = ~sum(collect_R2 < 0, 2);


collect_R2_val = collect_R2(non_negatives & non_nans, :); 

scatter_row_by_row(collect_R2, ylimits)

end
