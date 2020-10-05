
function [xVar] = bByb2HalvesXVariabilityPlot(bBybHalvesCell)
% bBybHalvesCell = stat.initHcell; 
    xVar = nan(size(bBybHalvesCell,1),4,2); % block-by-half 
    for b = 1:4 % the first four blocks
        for h = 1:2
            xVar(:,b,h) =  cell2mat(cellfun(@(a) nanstd(a(1,:))./sqrt(size(a,2)), bBybHalvesCell(:,b,h),'un',0)); 
        end 
    end
    plotIdvDataPointsTwoGroups(xVar);
end