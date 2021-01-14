function [jkvt,k] = evtKinematics( jkvt, sessionDur )

pullTrsIdx = cellfun(@(c)strcmpi(c,'sp'), {jkvt.trialType}); % successfull pull trials
k.jsVel   = zeros(1,sessionDur); % continuous js velocity
k.handVel = zeros(1,sessionDur); % continuous hand velocity

% For hand velocity, a realistic individual neuronal encoding might be reach vs. pull phases
% so first, just get hand trajectories on 1-ms resolution without binning
% not sure hTrjDstBase is the best measure to use - it might be. 
% bin next. Binning would do quite a bit of smoothing. 
% if not satisfactory, to further smoothing, as really the detailed
% velocity info is not the point here, it's really about going after
% distinguishing reach vs. pull phases (the directionality). 
% Of course I can try to assess encoding of velocity on each axis, but it's
% unlikely for individual neurons encoding that specifically. 

for t = 1:length(jkvt)
    % continuous hand velocity relative to the baseline point
    if ~isempty(jkvt(t).hTrjDstBase)
       handVel = diff([jkvt(t).hTrjDstBase(1) jkvt(t).hTrjDstBase]); % continuous hand kinetics
       hvT1 = jkvt(t).vFrameTime(1); % 1st time point of velocity trace
       hvTw = range(jkvt(t).vFrameTime)+1; % range
       intHandVel = interp1(jkvt(t).vFrameTime-hvT1+1, handVel, 1:hvTw,'linear'); % interpolate hand velocity
       % put zeros to points after the pullStop (they're due to joystick retraction)
       if ~isempty(jkvt(t).pullStops)
            intHandVel(min(jkvt(t).pullStops-hvT1+50,hvTw):end) = 0; % +50 give some room right after pullStop 
       end
       k.handVel(1,hvT1:hvT1+hvTw-1) = intHandVel./(10/1000); % convert to cm/s
    end
    
    if isfield(jkvt(t).movKins,'sgJsVel')
        k.jsVel(jkvt(t).trJsReady:jkvt(t).trJsReady+length(jkvt(t).movKins.sgJsVel)-1) = min(0,jkvt(t).movKins.sgJsVel);
    end
    
    % spot successful trials 
    if pullTrsIdx(t) && ~isempty(jkvt(t).movKins.pullStart) && ~isempty(jkvt(t).movKins.pullStop)
        % continuous js velocity take pulls (negative values) only
        %k.jsVel(jkvt(t).trJsReady:jkvt(t).trJsReady+length(jkvt(t).movKins.sgJsVel)-1) = min(0,jkvt(t).movKins.sgJsVel);
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