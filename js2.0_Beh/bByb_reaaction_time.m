function rt_blk = bByb_reaction_time(jkvt_blk, ss_blk)

rt_blk = nan(length(jkvt_blk), 1);


for jj = 1:length(jkvt_blk) 
    ttI = strcmpi(jkvt_blk(jj).trialType, 'sp'); 
    evtI = strcmpi(ss_blk(jj).evtAlign, 'rStart'); 

    if ttI && evtI
        if ~isempty(pullStart)
            rch_dur_blk(jj, 1) = max(0, pullStart - timeAlign); 
        end
        
        if ~isempty(pullStop)
            pull_dur_blk(jj, 1) = max(0, pullStop - pullStart); 
        end
    end
end
end
