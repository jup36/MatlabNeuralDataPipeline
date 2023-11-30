function flatCell = flattenNestCellArray(nestedCell)
flatCell = [];
for c = 1:length(nestedCell)
    flatCell = [flatCell; nestedCell{c}];
end
end