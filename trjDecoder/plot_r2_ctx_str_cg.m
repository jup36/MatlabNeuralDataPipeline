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