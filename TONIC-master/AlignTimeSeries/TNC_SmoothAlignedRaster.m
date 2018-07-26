function [smOut] = TNC_SmoothAlignedRaster(data,kernel,window)

tmp     = size(data,1);
smData  = zeros(size(data));

for i = 1:tmp

    smData(i,:) = conv(data(i,:),kernel,'same');    
    
end


smOut.image.psth           = sum(smData,1);
smOut.image.psthAVG        = mean(smData,1);
smOut.image.psthSEM        = std(smData,0,1) ./ sqrt(size(smData,1)-1);

% use the entire left window as an estimate of the mean
smOut.image.psthZ          = (smOut.image.psthAVG - mean(smOut.image.psthAVG(1:window(1,1)))) ./ std(smOut.image.psthAVG(1:window(1,1)));  
smOut.image.psthZe         = smOut.image.psthSEM ./ std(smOut.image.psthAVG(1:window(1,1)));