function [concatMat, numbUnit, numbTime, numbTrial] = concatUnitTimeBCell(unitTimeBCell) 
numbTrial = length(unitTimeBCell); 

sizeCell = cellfun(@size, unitTimeBCell, 'UniformOutput', false); 
sizeMat = cell2mat(sizeCell(:)); 
assert(length(unique(sizeMat(:, 1)))==1); % ensure that the number of units match
assert(length(unique(sizeMat(:, 2)))==1); % ensure that the number of time points match

numbUnit = unique(sizeMat(:, 1)); 
numbTime = unique(sizeMat(:, 2)); 

concatMat = full(cell2mat(cellfun(@(a) a', unitTimeBCell, 'UniformOutput', false)'));
end
