function [rch_dur_blk, pull_dur_blk] = bByb_reach_pull_duration(ss_blk)

rch_dur_blk = nan(length(ss_blk), 1);
pull_dur_blk = nan(length(ss_blk), 1);

for jj = 1:length(ss_blk) 
    ttI = strcmpi(ss_blk(jj).trialType, 'sp'); 
    evtI = strcmpi(ss_blk(jj).evtAlign, 'rStart'); 
    timeAlign = ss_blk(jj).timeAlign;
    pullStart = ss_blk(jj).tPullStart;
    pullStop = ss_blk(jj).tPullStop;
    
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
