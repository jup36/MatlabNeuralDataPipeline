function combineCell = combineCellsIntoCell(cellArr)
    % cellArr = rez.absRsAlignRchAngId; 
    combineCell = []; 
    for i = 1:length(cellArr)
        combineCell = [combineCell; cellArr{i}]; 
    end
end
