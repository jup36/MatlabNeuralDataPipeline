function [intMat] = intm(origMat, numbDataPoints)
% 1-d interpolation of a matrix as specified by the number of data points
    % numbDataPoints = 100; % 20ms*100 = 2000ms
    if numbDataPoints>size(origMat,2)
        x=1:size(origMat,2); 
        xq=linspace(1,size(origMat,2),numbDataPoints); 
        intMat = interp1(x,origMat',xq)'; 
    else
        intMat = origMat;
    end
end