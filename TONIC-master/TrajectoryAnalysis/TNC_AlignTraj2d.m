function [trajectories] = TNC_AlignTraj2d(x,y,v,stamps,window)

trajectories.x = zeros(numel(stamps),numel([-window(1):window(2)]));
trajectories.y = zeros(numel(stamps),numel([-window(1):window(2)]));
trajectories.v = zeros(numel(stamps),numel([-window(1):window(2)]));

for i=1:numel(stamps)
    
    if stamps(i) > window(1)
        if stamps(i) < numel(x)-window(2)
            trajectories.x(i,:) = x(stamps(i)-window(1):stamps(i)+window(2));
            trajectories.y(i,:) = y(stamps(i)-window(1):stamps(i)+window(2));
            trajectories.v(i,:) = v(stamps(i)-window(1):stamps(i)+window(2));
        end
    end
    
end