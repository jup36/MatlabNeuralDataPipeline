function [sR4c,sR5t] = sRateParse(sRateC)
% parse success rates by dividing each block into 4 chunks (chunk4c) or by dividing each block into sub-blocks of 5 trials 
    %sRateC = {rez.sRate}; 
    chunk4c = @(a) mat2cell(a(:), [fix(numel(a)/4) *[1 1 1], numel(a)-3*fix(numel(a)/4)], 1); % chuck the success logic array into 4 cells   
    chunk5t = @(a) mat2cell(a(:), [5*ones(1,fix(numel(a)/5)), numel(a)-5*fix(numel(a)/5)], 1); % chuck the success logic array into cells of 5 trials      
    
    for f = 1:length(sRateC) % file
        for b = 1:length(sRateC{f}) % block 
            si = sRateC{f}{b}; % this block's success logic
            if length(si)>=4
               sR4c{f,b} = cell2mat(cellfun(@(c) nansum(c)./length(c), chunk4c(si),'un',0));             
               sR5t{f,b} = cell2mat(cellfun(@(c) nansum(c)./length(c), chunk5t(si),'un',0));              
            else
               sR4c{f,b} = nan; 
               sR5t{f,b} = nan; 
            end
        end
    end
end