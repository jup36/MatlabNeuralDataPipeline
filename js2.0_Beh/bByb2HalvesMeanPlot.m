function [xMed] = bByb2HalvesMeanPlot(bBybHalvesCell)
% block-by-block variability of the initial hand position on X
% bBybHalvesCell = stat.initHcell; 
    xMed = nan(size(bBybHalvesCell,1),4,2); % block-by-half 
    for b = 1:4 % the first four blocks
        for h = 1:2
            xMed(:,b,h) =  cell2mat(cellfun(@(a) nanmean(a(1,:)), bBybHalvesCell(:,b,h),'un',0)); 
        end 
    end
    plotIdvDataPointsTwoGroups(xMed);
    hold off; 
end