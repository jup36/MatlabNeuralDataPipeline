function [map] = TNC_ConvertScatterToMatrix(data,dimensions,smthParams,scale)
% FUNCTION DETAILS: Use an X, Y vector pair of data and convert into an evenly sampled density matrix
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
%
% rescale the dimensions to fit on integer scaling if necessary (i.e. scale ~= 1)
    tmpData.x = data.x.*scale;
    tmpData.y = data.y.*scale;
    
% create a matrix to hold the full range of position data
    map.matrixForm = zeros( ceil(dimensions.y(1,2)-dimensions.y(1,1)+1) , ceil(dimensions.x(1,2)-dimensions.x(1,1)+1) );
    
% go through entire position data and accumulate samples per location
    for i=1:length(positions.x)
        map.matrixForm(round(tmpData.y(i)),round(tmpData.x(i))) = map.matrixForm(round(tmpData.y(i)),round(tmpData.x(i)))+1;
    end

% flexible scheme for smoothing the matrix form of the data with two-dimensional convolution with a gaussian kernel    
    if smthParams.sigma ~= 0
        map.smth.mu = smthParams.mu; % 16 is a reasonable default
        map.smth.sigma = smthParams.sigma;
        map.smth.gaussian = TNC_CreateGaussian(map.smth.mu,map.smth.sigma,map.smth.mu.*2,1);    
        map.smth.freq = conv2(map.matrixForm,map.smth.gaussian'*map.smth.gaussian,'same');
    else
        disp('Smoothing is not being used. Set sigma to a non-zero value to use smoothing of 2D maps');
        map.smth.mu = smthParams.mu; % 16 is a reaonable default
        map.smth.sigma = smthParams.sigma;
        map.smth.gaussian = 0;    
        map.smth.freq = spikeMap.freq;
    end
