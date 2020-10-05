function [chk4c,chk5t,chk2c] = blockByBlockChunkMean(bBybCell)
%This function chuncks each block of a block-by-block cell into n subblocks or subblocks of m trials and gets the mean of the each subblock   
    chunk4c = @(a) mat2cell(a(:), [fix(numel(a)/4) *[1 1 1], numel(a)-3*fix(numel(a)/4)], 1); % chuck the success logic array into 4 cells   
    chunk5t = @(a) mat2cell(a(:), [5*ones(1,fix(numel(a)/5)), numel(a)-5*fix(numel(a)/5)], 1); % chuck the success logic array into cells of 5 trials     
    chunk2c = @(a) mat2cell(a(:), [fix(numel(a)/2), numel(a)-fix(numel(a)/2)], 1);
    for f = 1:length(bBybCell) % file
        for b = 1:length(bBybCell{f}) % block 
            si = bBybCell{f}{b}; % this block's success logic
            if iscell(si) 
                si = cell2mat(si); 
            end
            if sum(isinf(si))>0
                si(isinf(si))=nan;
            end 
            if length(si)>=4
               chk4c{f,b} = cell2mat(cellfun(@nanmean, chunk4c(si),'un',0));             
               chk5t{f,b} = cell2mat(cellfun(@nanmean, chunk5t(si),'un',0));
               chk2c{f,b} = cell2mat(cellfun(@nanmean, chunk2c(si),'un',0));                
            else
               chk4c{f,b} = nan; 
               chk5t{f,b} = nan; 
               chk2c{f,b} = nan;  
            end
        end
    end
end
