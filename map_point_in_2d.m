function minI = map_point_in_2d(point)
ref = -2:.2:2; 
[~, minI] = min(abs(point - ref));

end