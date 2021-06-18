function c = emptyCell2nan(c)      %named function
  c(cellfun(@isempty, c)) = {nan};
end
