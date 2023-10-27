function [lickCnt, lickInt, lickLambda, lickIntWhole] = examineLickC(lickC, minMaxTime)

lickCnt = nanmean(cell2mat(cellfun(@(a) length(a)/abs(diff(minMaxTime)), lickC, 'UniformOutput', false))); 
lickInt = nanmean(cell2mat(cellfun(@(a) diff(sort([minMaxTime(1); a])), lickC, 'UniformOutput', false)'));
lickIntWhole = cell2mat(cellfun(@(a) diff(sort([minMaxTime(1); a])), lickC, 'UniformOutput', false)'); 
lickLambda = 1/lickInt;

end