function [jkvt,k] = evtKinematics( jkvt, sessionDur )

pullTrsIdx = cellfun(@(c)strcmpi(c,'sp'), {jkvt.trialType}); % successfull pull trials
k.jsVel   = zeros(1,sessionDur); % continuous js velocity
k.handVel = zeros(1,sessionDur); % continuous hand velocity

for t = 1:length(jkvt)
    % continuous hand velocity relative to the baseline point
    if ~isempty(jkvt(t).hTrjDstBase)
       handVel = diff([jkvt(t).hTrjDstBase(1) jkvt(t).hTrjDstBase])./(1/1000); % continuous hand kinetics
       intHandVel = interp1(jkvt(t).vFrameTime-jkvt(t).vFrameTime(1)+1, handVel, 1:range(jkvt(t).vFrameTime)+1,'linear'); % interpolate hand velocity
       k.handVel(1,jkvt(t).vFrameTime(1):jkvt(t).vFrameTime(1)+length(intHandVel)-1) = intHandVel; 
    end
    
    % spot successful trials 
    if pullTrsIdx(t) && ~isempty(jkvt(t).movKins.pullStart) && ~isempty(jkvt(t).movKins.pullStop)
        % continuous js velocity take pulls (negative values) only
        k.jsVel(jkvt(t).trJsReady:jkvt(t).trJsReady+length(jkvt(t).movKins.sgJsVel)-1) = min(0,jkvt(t).movKins.sgJsVel);
        % detect pull start/stop of successful trials
        jkvt(t).pullStarts = jkvt(t).trJsReady + jkvt(t).movKins.pullStart; 
        jkvt(t).pullStops  = jkvt(t).trJsReady + jkvt(t).movKins.pullStop;         
        % detect reach start (hand lift) of each successful trial 
        tmphTrjRstartT = jkvt(t).vFrameTime(jkvt(t).hTrjRstart); 
        jkvt(t).rStartToPull = tmphTrjRstartT(find(tmphTrjRstartT<jkvt(t).pullStarts,1,'last'));
        % detect reach stop 
        tmphTrjRstopT = jkvt(t).vFrameTime(jkvt(t).hTrjRstop); 
        jkvt(t).rStopToPull = tmphTrjRstopT(find(tmphTrjRstopT>jkvt(t).pullStarts,1,'last'));        
    end
end

end