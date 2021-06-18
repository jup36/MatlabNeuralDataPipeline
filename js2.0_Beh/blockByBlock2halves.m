function [bBybHalvesC,bBybHalvesM,bBybHalvesSem] = blockByBlock2halves(bBybCell,bBybRwdCell)
%This function takes an array of cells, each element of which are block-by-block organized, and partition each block's valid trials data into two groups 
% bBybCell = {rez.hXYtort}; 
% bBybRwdCell = {rez.rwdTrI}; 
for f = 1:length(bBybCell)
  for b = 1:length(bBybCell{f})
      %c = bBybCell{f}{b}; 
      c = cellfun(@(a) a(:,end), bBybCell{f}{b}, 'un', 0); 
      r = bBybRwdCell{f}{b}; 
      rValTrs = find(r&~cell2mat(cellfun(@(a) sum(isinf(a))>0|sum(isnan(a))>0,c,'un',0))); % index of rewarded and valid (non-NaN) trials 
      % should I also do valid trials halves regardless of reward?
      if length(rValTrs)>=4 % if there're more than 4 rewarded trials in the block
          %hf1 = rValTrs(1:fix(sum(r)/2)); % 1st half trials 
          %hf2 = rValTrs(fix(sum(r)/2)+1:end); % 2nd half trials
          hf1 = rValTrs(1:fix(length(rValTrs)/2)); % 1st half trials 
          hf2 = rValTrs(fix(length(rValTrs)/2)+1:end); % 2nd half trials
          bBybHalvesC{f,b,1} = cell2mat(c(hf1));  
          bBybHalvesC{f,b,2} = cell2mat(c(hf2));
          bBybHalvesM{f,b,1} = nanmedian(bBybHalvesC{f,b,1}'); 
          [~,~,bBybHalvesSem{f,b,1}] = meanstdsem(bBybHalvesC{f,b,1}'); 
          bBybHalvesM{f,b,2} = nanmedian(bBybHalvesC{f,b,2}'); 
          [~,~,bBybHalvesSem{f,b,2}] = meanstdsem(bBybHalvesC{f,b,2}'); 
      end
  end   
end
end
