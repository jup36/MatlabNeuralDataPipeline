function refTsI = check_timestamp_overlap_adjacentFrames(refTs, evtTs, nAdjacent)
% Generate example data (replace with your actual data)
%refTs = faceCam_pt;
%evtTs = xStampSec_lick;

% Find the minimum and maximum timestamps of the first set
minRefTs = min(refTs);
maxRefTs = max(refTs);

refTsI = zeros(length(refTs), 1); 

% Select second timestamps that fall between the minimum and maximum of the first set
selectEvtTs = evtTs(evtTs(:, 1) >= minRefTs & evtTs(:, 1) <= maxRefTs, 1);

if ~isempty(selectEvtTs)
    for ts = 1:length(selectEvtTs)
        [~, minI] = min(abs(refTs - selectEvtTs(ts, 1)));
        refTsI(max(1, minI-nAdjacent):min(minI+nAdjacent, length(refTsI))) = 1;
    end
end

end