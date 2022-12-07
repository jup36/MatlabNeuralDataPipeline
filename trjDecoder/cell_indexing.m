function c_indexed = cell_indexing(c_to_index, index)
c_indexed = cell(size(c_to_index));
val_c_I = cell2mat(cellfun(@(a) ~isempty(a), c_to_index, 'un', 0));
val_c_indexed = cellfun(@(a) a(index, :), c_to_index(val_c_I), 'un', 0);
c_indexed(val_c_I) = deal(val_c_indexed);
end