function [excludeLogic] = excludeSegments(dat,pointsToExclude,oneWindow)
% this function takes a whole segment and excludes
%subsegment-windows specified by the points and window.
% dat = X;
% pointsToExclude = testPtsToExclude_bin;
% oneWindow = 2000/p.DT;

if size(dat,1)<size(dat,2)
    dat = dat';
end

windowsToExclude = arrayfun(@(a) a-oneWindow:a+oneWindow, pointsToExclude, 'un', 0);
excludeLogicC = cellfun(@(a) ismember(1:size(dat,1),a), windowsToExclude, 'un', 0)';
excludeLogic = sum(cell2mat(excludeLogicC))>0;

if size(excludeLogic,1)<size(excludeLogic,2)
    excludeLogic = excludeLogic';
end

end
