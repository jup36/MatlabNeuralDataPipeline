function [evtAlignedDffTs, evtAlignedTs] = alignToEvent3D(dffTs, eventTimeToAlign, frameT, timeWin)
% NOTE: Modified the function to return outcomes even if the timeWin goes
% out of bounds of the frame time! (2/4/2025, Junchol Park)

win = timeWin + eventTimeToAlign;

frameI = frameT >= min(win) & frameT <= max(win);
evtAlignedDffTs = dffTs(:, :, frameI);
evtAlignedTs = frameT(frameI);

end