function [hitLickLatencyC] = hitLickRastersAuditoryGng(filePath)

[~, header] = fileparts(filePath);

% get behavioral data further analyzed
fileBehParseGng = GrabFiles_sort_trials('tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehParseGng)
    tbytDat = parseAuditoryGngTrials(tbytDat);
    save(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')
else
    load(fileBehParseGng{1}, 'tbytDat')
end

%% indices for trial-type identification
waterI = cellfun(@(a) ~isempty(a), {tbytDat.water}); 
lickI = cellfun(@(a) ~isempty(a), {tbytDat.Lick}); 
airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff}); 
goI = [tbytDat.rewardTrI]'==1; 
nogoI = [tbytDat.punishTrI]'==1; 
hitI = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})'; 
missI = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})'; 
faI = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})'; 
crI = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})'; 

%
hitLickLatencyC = cellfun(@(a, b) a-b, {tbytDat(hitI).Lick}, {tbytDat(hitI).evtOn}, 'UniformOutput', false); 


end