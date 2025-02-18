function [dffOut, tF] = alignToEventInterp3D(dffTs, eventTimeToAlign, frameT, timeWin, step)
% NOTE: Modified the function to return outcomes even if the timeWin goes
% out of bounds of the frame time! (2/4/2025, Junchol Park)
% e.g., step = 0.05; % 50ms

alignedFrameT = frameT + eventTimeToAlign;  
tF = timeWin(1):step:timeWin(end); 

% Define the valid range based on frameTOnTf
validTfIdx = tF >= min(alignedFrameT) & tF <= max(alignedFrameT);
tFval = tF(validTfIdx); 

% Get dimensions
[rows, cols, ~] = size(dffTs);

% Initialize output with NaNs 
dffOut = NaN(size(dffTs,1),size(dffTs,2),length(tF)); 

% Loop through each pixel
for i = 1:rows
    for j = 1:cols
        % Extract time series for this pixel
        pixelSeries = squeeze(dffTs(i, j, :));
        
        % Apply strict NaN criterion: Only interpolate if ALL time points are valid
        if all(~isnan(pixelSeries))
            dffOut(i, j, validTfIdx) = interp1(alignedFrameT, pixelSeries, tFval, 'pchip');
        end
    end
end

end