function outcell = plot_r2_ctx_str_cg(r2_cell, ylimits, col)

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
outcell = cell(size(collect_R2_C, 1), size(collect_R2_C, 2)); 
non_nans = cell2mat(cellfun(@(a) ~isempty(sum(a, 2)), collect_R2_C, 'un', 0)); 
[outcell{~non_nans}] = deal(NaN); 

collect_R2_val = cellfun(@(a) a(col), collect_R2_C(non_nans), 'un', 0); 

outcell(non_nans) = collect_R2_val; 
outcell(cell2mat(cellfun(@(a) a<0, outcell, 'un', 0))) = {NaN}; % take negative values as NaN

scatter_row_by_row_cell(outcell, ylimits)

end