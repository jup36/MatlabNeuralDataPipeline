function p = parse_input_jsVideo(filePath, nTr, vargs)
% parse input, and extract name-value pairs
default_frameRate = 250; % the default frame rate of videos
default_slowPlay = 1; % fold to be slowed down for playback
default_frTimeRange = []; % the frame time range to be included in the movie, if [] - include all. To select frames, e.g. [-1000 500], the frames within the 1500 ms range relative to the event will be selected.

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'nTr')
addParameter(p,'frameRate',default_frameRate)
addParameter(p,'slowPlay',default_slowPlay)
addParameter(p,'frTimeRange',default_frTimeRange)

parse(p,filePath,nTr,vargs{:})
end