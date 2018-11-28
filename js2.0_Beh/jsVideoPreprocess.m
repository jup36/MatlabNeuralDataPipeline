function jsVideoPreprocess(filePath, varargin)

filePath = 'C:\Users\parkj\Documents\DATA\VideoAPT\WR25_111618';
cd(filePath)

p = parse_input_jsVideo(filePath, varargin); 
% p = parse_input_jsVideo(filePath, {}); 

close all; clc

vid1 = VideoReader('cam0_cam_0_2018_11_16_15_40_37.avi');
vid2 = VideoReader('cam1_cam_1_2018_11_16_15_40_37.avi');

%videoPlayer = vision.VideoPlayer;

% new video
outputVideo = VideoWriter('newvideo.avi');
outputVideo.FrameRate = vid1.FrameRate;
open(outputVideo);

while hasFrame(vid1) && hasFrame(vid2)
    img1 = readFrame(vid1);
    img2 = readFrame(vid2);

    imgt = horzcat(img1, img2);

    % play video
    %step(videoPlayer, imgt);

    % record new video
    writeVideo(outputVideo, imgt);
end

%release(videoPlayer);
close(outputVideo);

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_jsVideo(filePath, vargs)
        % parse input, and extract name-value pairs
        default_frameRate = 250; % the default frame rate of videos
        default_slowPlay = 5; % fold to be slowed down for playback
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'frameRate',default_frameRate)
        addParameter(p,'slowPlay',default_slowPlay)
        
        parse(p,filePath,vargs{:})    
    end
end