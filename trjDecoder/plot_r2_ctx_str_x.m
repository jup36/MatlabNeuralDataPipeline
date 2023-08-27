function collect_R2_val = plot_r2_ctx_str_x(r2_cell, ylimits)

ctxI = cell2mat(cellfun(@(a) strcmp(a, 'ctx'), r2_cell(1,:), 'un', 0)); 
strI = cell2mat(cellfun(@(a) strcmp(a, 'str'), r2_cell(1,:), 'un', 0)); 
ctx_strI = cell2mat(cellfun(@(a) strcmp(a, 'ctx_str'), r2_cell(1,:), 'un', 0)); 

ctx_R2 = r2_cell(2:end, ctxI); 
str_R2 = r2_cell(2:end, strI); 
ctx_str_R2 = r2_cell(2:end, ctx_strI); 

collect_R2_C = [ctx_R2, str_R2, ctx_str_R2]; 
empties = cell2mat(cellfun(@isempty, collect_R2_C, 'un', 0));

[collect_R2_C{empties}] = deal(NaN); 

non_nans = cell2mat(cellfun(@(a) ~isnan(sum(a, 2)), collect_R2_C, 'un', 0)); 
non_negatives = cell2mat(cellfun(@(a) sum(a < 0)==0, collect_R2_C, 'un', 0)); 
collect_R2 = collect_R2_C(prod(non_negatives & non_nans, 2)==1, :); 
collect_R2_val = cell2mat(cellfun(@(a) a(1), collect_R2, 'un', 0)); 
scatter_row_by_row(collect_R2_val, ylimits)

%non_nans = ~isnan(sum(collect_R2, 2)); 
%non_negatives = ~sum(collect_R2 < 0, 2);
%collect_R2_val = collect_R2(non_negatives & non_nans, :); 
%collect_R2 = cellfun(@(a) a(:), collect_R2_C, 'un', 0); %cell2mat(cellfun(@(a) a(:), collect_R2_C, 'un', 0)); 
end
