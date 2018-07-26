function [geometry] = getHH2AcuteGeom
    vi32 = 1:32;
    geometry = zeros(64, 2);
    geometry(vi32,2) = 25*(vi32-1);
    geometry(vi32+32,1) = 250;
    geometry(vi32+32,2) = geometry(vi32,2);
end