function minI = map_point_X(point)
ref = -2:.2:2; 
[~, minI] = min(abs(point - ref));
end