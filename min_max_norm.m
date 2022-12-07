function x_n = min_max_norm(x, min_sub)
x_m = nanmean(x, 1); 

if min_sub == true 
    x_m_min = min(x_m); 
    x_m = x_m - x_m_min; 
end

x_m_max = max(x_m);
x_n = x_m/x_m_max; 
end