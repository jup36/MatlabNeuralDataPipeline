function [spd_blk, spd_blk_x, spd_blk_y, spd_blk_z] = bByb_maxRchSpeed(ss_blk)

spd_blk = nan(length(ss_blk), 1);

for jj = 1:length(ss_blk) 
    ttI = strcmpi(ss_blk(jj).trialType, 'sp'); 
    evtI = strcmpi(ss_blk(jj).evtAlign, 'rStart'); 

    if ttI && evtI
        if ~isempty(ss_blk(jj).maxRchSpeed) && ~isempty(ss_blk(jj).maxRchSpeedXYZ)
            spd_blk(jj, 1) = ss_blk(jj).maxRchSpeed;
            spd_blk_x(jj, 1) = ss_blk(jj).maxRchSpeedXYZ(1);
            spd_blk_y(jj, 1) = ss_blk(jj).maxRchSpeedXYZ(2);
            spd_blk_z(jj, 1) = ss_blk(jj).maxRchSpeedXYZ(3);
    
        end

    end
end
end
