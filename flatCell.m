function flatC = flatCell(nestedCell)
flatC = [];
for c = 1:length(nestedCell)
    if ~isempty(nestedCell{c})
        flatC = [flatC; nestedCell{c}];

    end
end
end