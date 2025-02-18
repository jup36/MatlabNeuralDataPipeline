function index = taskEventIndices(tbytDat) 
%This utility function generates various task event indices from the tbytDat structure. 

index.waterI = cellfun(@(a) ~isempty(a), {tbytDat.water}); 
index.lickI = cellfun(@(a) ~isempty(a), {tbytDat.Lick}); 
index.airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff}); 
index.goI = [tbytDat.rewardTrI]'==1; 
index.nogoI = [tbytDat.punishTrI]'==1; 
index.hitI = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})'; 
index.missI = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})'; 
index.faI = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})'; 
index.crI = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})'; 
end