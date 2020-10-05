function [yVar] = bByb2HalvesYVariabilityPlot(bBybHalvesCell)
% bBybHalvesCell = stat.initHcell; 
    yVar = nan(size(bBybHalvesCell,1),4,2); % block-by-half 
    for b = 1:4 % the first four blocks
        for h = 1:2
            yVar(:,b,h) =  cell2mat(cellfun(@(a) nanstd(a(2,:))./sqrt(size(a,2)), bBybHalvesCell(:,b,h),'un',0)); 
        end 
    end
    plotIdvDataPointsTwoGroups(yVar);
end
