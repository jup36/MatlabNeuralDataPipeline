function flatC = flatCell(nestedCell)
flatC = [];
for c = 1:length(nestedCell)
    flatC = [flatC; nestedCell{c}];
end
end