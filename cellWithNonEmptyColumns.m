function [nonEmptyCell] = cellWithNonEmptyColumns(cellArray)

numCol = size(cellArray, 2);  

valRowI = sum(cell2mat(cellfun(@(a) ~isempty(a), cellArray, 'UniformOutput', false)), 2)==numCol; 

nonEmptyCell = cellArray(valRowI, :); 


end