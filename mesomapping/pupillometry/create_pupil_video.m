function create_pupil_video(file_path, video_name, mID, date, likelihoodcut, frameRate)

%% collect info
file_path = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/cj942/m942_070623';
video_name = 'behvid_13_22_52_942_2_cropped.mp4';
mID = '942';
date = '070623';
likelihoodcut = 0.5; 
frameRate = 400; 

%% load stimopts
filePath_stimopts = findFileWithString(file_path, 'stimInfo');
load(filePath_stimopts)

%% load nidq data
filePath_nidq = findFolderWithString(file_path, '_g');
if isempty(filePath_nidq)
    filePath_nidq = uigetdir(file_path, 'Select a folder');
end

timestamp_behav_events(filePath_nidq, 'faceCam', 'lick')

load(fullfile(filePath_nidq, 'evtInS.mat'), 'evtInS')

%% load trial-by-trial behav data
behFolder = dir('/Volumes/buschman/Users/Caroline/NADA_dynamics/data');
tbytI = cell2mat(cellfun(@(a) contains(a, mID) && contains(a, date), {behFolder.name}, 'UniformOutput', false));
load(fullfile(behFolder(tbytI).folder, behFolder(tbytI).name), 'data', 'xStampSec_airpuff', 'xStampSec_lick', 'xStampSec_water');
tbytDat = data;
clearvars behFolder data

%% get video file path and csv file path
filePath_vid = findFolderWithString(file_path, '_vid_cropped');
if isempty(filePath_vid)
    warning("The path for cropped videos doesn't exist!")
else
    filePath_vidFile = fullfile(filePath_vid, video_name);
end
[~, video_name_base, ~] = fileparts(video_name);

filePath_raw_vid = findFolderWithString(file_path, '_vid');
if isempty(filePath_raw_vid)
    warning("The path for raw videos doesn't exist!")
else
    filePath_rawVidFile = fullfile(filePath_raw_vid, [video_name_base(1:end-8), '.mp4']);
end

% read csv file with coordinates
filePath_csv = findCSVFileWithString(filePath_vid, video_name_base);
if ~isempty(filePath_csv)
    fprintf('Found the csv file for cooridnates!\n');
end

tableDlc = readDlcCsv(filePath_csv);

%% get behavioral events relative to frames
% get the number of trial
vReaderP = VideoReader(filePath_vidFile);
vReaderF = VideoReader(filePath_rawVidFile);

matches = regexp(video_name, '_(\d+)_cropped', 'tokens');
numTr = str2double(matches{1}{1});

faceCam_pt = evtInS.faceCam(evtInS.faceCam(:, 2)==numTr, 1);

if length(faceCam_pt) > vReaderP.NumFrames
    faceCam_pt = faceCam_pt(1:vReaderP.NumFrames);
end

%% map temporal events to camera pulses
frameLickI = check_timestamp_overlap(faceCam_pt, xStampSec_lick(:, 1));
frameWaterI = check_timestamp_overlap(faceCam_pt, xStampSec_water(:, 1));
frameAirpuffI = check_timestamp_overlap(faceCam_pt, xStampSec_airpuff(:, 1));

% timestamp each frame relative to the stim onset
frameStimOnI = check_timestamp_overlap(faceCam_pt, tbytDat(numTr).start_stim);
frameStimOffI = check_timestamp_overlap(faceCam_pt, tbytDat(numTr).end_stim);
frameStimI = zeros(length(frameStimOnI), 1);
frameStimI(find(frameStimOnI, 1):find(frameStimOffI, 1), 1) = 1;
reFaceCam_pt = round(faceCam_pt-tbytDat(numTr).start_stim, 2);

%% generate the labeled video
% first read and label frames and save
assert(vReaderP.NumFrames == size(tableDlc, 1));

dotColors = {[120, 95, 170], [81, 121, 235], [132, 171, 181], [149, 244, 207], ...
    [192, 249, 171], [231, 203, 128], [230, 131, 75], [206, 57, 32]}; % from the left clockwise
dotColors = cellfun(@(a) a./255, dotColors, 'un', 0);

filePath_png = fullfile(file_path, strcat(video_name_base, '_png'));

if exist(filePath_png, 'dir') == 0
    mkdir(filePath_png);
end

for ff = 1:vReaderP.NumFrames
    frameP = readFrame(vReaderP);
    frameF = readFrame(vReaderF);
    frameF = frameF(:, 1:360, :); % cut the facial half of the frameF

    % enlarge the pupil frame to match the dimension (height)
    frameP_lg = imresize(frameP, [size(frameF, 1), NaN]);

    % calculate the resizing factor and transform the dot coordinates
    resizeF = size(frameF, 1)/size(frameP, 1);

    dots_raw = get_pupil_dots(tableDlc(ff, :), likelihoodcut);
    dots_raw = cell2mat(dots_raw);
    dots = dots_raw * resizeF;

    % draw the frame
    combined = [frameP_lg, frameF];
    imshow(combined); hold on;

    % draw dots
    for dt = 1:size(dots, 1)
        if ~isnan(dots(dt, 1))
            % Add dots to the frame with different colors and transparency
            scatter(dots(dt, 1), dots(dt, 2), 150, 'MarkerFaceColor', dotColors{dt}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5);
        end
    end

    % texts
    textPos = {[20, 410], [190, 410], [380, 410]};

    textStr = cell(1, 3);
    textColor = 'white';
    fontSize = 14;

    % timer
    if ~isempty(reFaceCam_pt(ff))
        textStr{1, 1} = sprintf('%.2fs', reFaceCam_pt(ff));
    end

    % lick tracker
    if frameLickI(ff)==1 % if there was a lick
        textStr{1, 2} = 'Lick';
    end

    % water tracker
    if frameWaterI(ff)==1 % if there was a drop
        textStr{1, 3} = 'Water';
    end

    % airpuff tracker
    if frameAirpuffI(ff)==1 % if there was an airpuff
        textStr{1, 3} = 'Airpuff';
    end

    for i = 1:numel(textPos)
        str = textStr{1, i};
        if ~isempty(str)
            pos = textPos{1, i};

            % Create a semi-transparent rectangle
            %rectanglePosition = [pos(1) - 5, pos(2)-60, 105, 47]; % x y w h
            %rectangleColor = [0, 0, 1, 0.3];  % Transparent black color (R, G, B, Alpha)

            % Add the rectangle
            %rectangle('Position', rectanglePosition, 'EdgeColor', 'none', 'FaceColor', rectangleColor);

            % Add the text
            text(pos(1), pos(2), str, 'Color', textColor, 'FontSize', fontSize);
        end
    end
    hold off;
    % Capture the current figure with dots
    frameLabled = getframe(gcf);
    frameLabled = cropborder(frameLabled.cdata,[NaN NaN NaN NaN],'threshold',0.0001); % from MIMT
    %frameLabled = frameLabled.cdata;

    % save the figure with dots
    png_name = sprintf('concat_frame_%d.png', ff);
    imwrite(frameLabled, fullfile(filePath_png, png_name));

    fprintf('Frame #%d is labeled and saved.\n', ff);
end


labeled_videos = findAsManyFilesWithString(filePath_vid, strcat(video_name_base, '_', 'labeled')); 

labeledVideo = VideoWriter(fullfile(filePath_vid, [strcat(video_name_base, '_', 'labeled'), sprintf('_%d', length(labeled_videos)), sprintf('_FR%d', frameRate)]), 'MPEG-4');
labeledVideo.FrameRate = frameRate;
playspeed = labeledVideo.FrameRate/200; 
open(labeledVideo);

for ff = 1:length(pngFiles)
    % Load PNG image
    imagePath = fullfile(filePath_png, sprintf('concat_frame_%d.png', ff));
    img = imread(imagePath); 
    
    % Add text to the image
    textStr = sprintf('x%.1f', playspeed);
    imgWithText = insertText(img, [20, 360], textStr, 'FontSize', 18, 'TextColor', [226, 205, 132]./255, 'BoxOpacity', 0.4);

    frame = im2frame(imgWithText); 

    % write the labeled frame to the output video
    writeVideo(labeledVideo, frame);
    close; 
    fprintf('Frame #%d is loaded and written.\n', ff);

end

% close the output video
close(labeledVideo);

disp('Video creation complete.');

end






