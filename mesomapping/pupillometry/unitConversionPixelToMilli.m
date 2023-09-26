function pixelToMm = unitConversionPixelToMilli(videoFilePath)

% Create a VideoReader object to read the video
videoReader = VideoReader(videoFilePath);

% Read the first frame
firstFrame = readFrame(videoReader);

[point1, point2] = getUserDrawnLineWithFeedback(firstFrame);
knownLength_pixels = abs(point2(2)-point1(2));
knownLength_mm = 2; % 2 mm (the metal platform to which the headplate is fixed)

% Calculate conversion factor for length
pixelToMm = knownLength_mm / knownLength_pixels;

save(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj', 'pixelToMm'), 'pixelToMm');
end


