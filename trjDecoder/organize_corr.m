function corrRez_collect = organize_corr(corrRez, corrRez_collect, row_I, cell_I, saveName)

fields = fieldnames(corrRez); 

corrRez_collect{row_I, 1} = saveName;
for j = 1:length(fields)
    field_name = fields{j}; 
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), corrRez_collect(1,:), 'un', 0)); 
    if isstruct(corrRez.(field_name){1,cell_I})
        corr_xyz = diag(corrRez.(field_name){1,cell_I}.all_sm)'; 
    else
        corr_xyz = []; 
    end
    corrRez_collect{row_I, col_I} = corr_xyz; 
end

end