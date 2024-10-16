function cmosExposureEventsC = organize_cmos_exposure_pulses(cmosExposureEvents, frameShortageAllowance)
% cmosExposureEvents is assumed to be two columns:
%   1) time in sec each pulse
%   2) trial indices

cmosExposureEventsC = cell(length(unique(cmosExposureEvents(:, 2))), 1);
for t = 1:length(unique(cmosExposureEvents(:, 2)))
    cmosExposureEventsC{t} = cmosExposureEvents(cmosExposureEvents(:, 2)==t, 1);
end

meanFrames = mean(cell2mat(cellfun(@length, cmosExposureEventsC, 'un', 0))); 
enoughFrameLogic = cell2mat(cellfun(@(a) length(a)>=meanFrames-frameShortageAllowance, cmosExposureEventsC, 'un', 0));

cmosExposureEventsC = cmosExposureEventsC(enoughFrameLogic); 

end