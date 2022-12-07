function r2Rez_collect = organize_r2_pos(r2Rez, r2Rez_collect, row_I, cell_I)

fields = fieldnames(r2Rez); 

for j = 1:length(fields)
    field_name = fields{j}; 
    r2_xyz_all = [r2Rez.(field_name){1,cell_I}.all, r2Rez.(field_name){1,cell_I}.overall]; 
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), r2Rez_collect(1,:), 'un', 0)); 
    r2Rez_collect{row_I, col_I} = r2_xyz_all; 
end

end