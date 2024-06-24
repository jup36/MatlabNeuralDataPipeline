    function trialI = trialTypesFromTbytDatParseGng(tbytDat)
        trialI.hitI  = [tbytDat.rewardTrI] & cellfun(@(a) ~isempty(a), {tbytDat.hitLicks}, 'UniformOutput', true); 
        trialI.missI = [tbytDat.rewardTrI] & cellfun(@(a) isempty(a), {tbytDat.hitLicks}, 'UniformOutput', true); 
        trialI.faI = [tbytDat.punishTrI] & cellfun(@(a) ~isempty(a), {tbytDat.faLicks}, 'UniformOutput', true); 
        trialI.crI = [tbytDat.punishTrI] & cellfun(@(a) isempty(a), {tbytDat.faLicks}, 'UniformOutput', true); 
    end