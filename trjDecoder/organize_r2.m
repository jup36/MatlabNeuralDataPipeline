function r2Rez_collect = organize_r2(r2Rez, r2Rez_collect, row_I, cell_I, saveName)

fields = fieldnames(r2Rez);

r2Rez_collect{row_I, 1} = saveName;
for j = 1:length(fields)
    field_name = fields{j};
    
    col_I = cell2mat(cellfun(@(a) strcmp(a, field_name), r2Rez_collect(1,:), 'un', 0));
    if isstruct(r2Rez.(field_name){1,cell_I})
        r2_xyz_all = [r2Rez.(field_name){1,cell_I}.all, r2Rez.(field_name){1,cell_I}.overall];
    else
        r2_xyz_all = [];
    end
    r2Rez_collect{row_I, col_I} = r2_xyz_all;
end

end