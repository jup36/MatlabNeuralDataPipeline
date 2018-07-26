function X_row = get_row_means(block_struct)
% for blockproc, compute means along rows of distance matrix

X_row = mean(block_struct.data,2);
X_row = repmat(X_row, 1, size(block_struct.data,2));
end

