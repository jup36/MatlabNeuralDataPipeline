function rt_blk = bByb_reaction_time(jkvt_blk, ss_blk)

rt_blk = nan(length(jkvt_blk), 1);


for jj = 1:length(jkvt_blk) 
    ttI = strcmpi(jkvt_blk(jj).trialType, 'sp'); 
    evtI = strcmpi(ss_blk(jj).evtAlign, 'rStart'); 

    if ttI && evtI
       rt_blk(jj, 1) = min(jkvt_blk(jj).rt, 10000); 
    end
end
end
