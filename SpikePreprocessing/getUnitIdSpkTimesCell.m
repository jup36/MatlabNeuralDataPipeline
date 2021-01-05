function [unitIdSTC] = getUnitIdSpkTimesCell(spkTimesCellIn, sIn, unitIdS)
%unit Ids don't match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
% this function identifies the unit from 'S' in 'spkTimesCell' and outputs
% it as 'unitidSTC'

geomI = cell2mat(cellfun(@(a) isequal(a,sIn.geometry{unitIdS}), spkTimesCellIn(4,:), 'un', 0)); 
wfI = cell2mat(cellfun(@(a) isequal(a,sIn.meanWF{unitIdS}), spkTimesCellIn(6,:), 'un', 0)); 

unitIdSTC = find(geomI & wfI); 

end