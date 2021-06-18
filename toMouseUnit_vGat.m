function [merged] = toMouseUnit_vGat(fileList,fileToMerge)
    %fileToMerge=rwdPosCMean;  
    %fileList = eachFile;     
    orgList = cellfun(@(a) a(1:6), fileList, 'un', 0); 
    uniquelist = unique(cellfun(@(a) a(1:6), fileList, 'un', 0)); 
    merged = nan(length(uniquelist),size(fileToMerge,2));   
    for f = 1:length(uniquelist)
        merged(f,:) = nanmean(fileToMerge(cell2mat(cellfun(@(a) strcmpi(a,uniquelist{f}), orgList, 'un', 0)), :),1);           
    end
end