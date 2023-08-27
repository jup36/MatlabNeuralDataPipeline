
filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/cj935/m935_062823/m935_062823_vid';
fileName = 'behvid_18_45_41_m935_20.mp4'; 

% video file save dir
tbytVideoPath = fullfile(filePath,'tbytVideos');
if ~isfolder(fullfile(filePath,'tbytVideos'))
    mkdir(fullfile(filePath,'tbytVideos'))
end

vid = VideoReader(fullfile(filePath, fileName));
