function corrRez_collect = organize_corr_pos(corrRez, corrRez_collect, row_I, cell_I)

fields = fieldnames(corrRez); 

for j = 1:length(fields)
    field_name = fields{j}; 
    corr_xyz = diag(corrRez.(field_name){1,cell_I}.all_sm)'; 
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), corrRez_collect(1,:), 'un', 0)); 
    corrRez_collect{row_I, col_I} = corr_xyz; 
end

end