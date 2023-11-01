function collect_R2_val = plot_r2_ctx_str_cg_pull_Z(r2_cell, ylimits)

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

[collect_R2_C{empties}] = deal([NaN, NaN, NaN, NaN]); 

non_nans = cell2mat(cellfun(@(a) ~isnan(sum(a, 2)), collect_R2_C, 'un', 0)); 
collect_R2 = collect_R2_C(sum(non_nans, 2)>=1, :); 

collect_R2_val = cell2mat(cellfun(@(a) a(3), collect_R2, 'un', 0)); 
val_I = prod(~isnan(collect_R2_val), 2); 
collect_R2_val_ = collect_R2_val(val_I==1, :); 
collect_R2_val_ = [collect_R2_val_; nansum(collect_R2_val(val_I~=1, :))];

scatter_row_by_row(collect_R2_val_, ylimits)


end