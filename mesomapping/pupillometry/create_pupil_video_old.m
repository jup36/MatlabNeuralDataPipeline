function create_pupil_video(file_path, video_name, mID, date, likelihoodcut)

%% collect info
% file_path = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/cj942/m942_070623';
% video_name = 'behvid_13_22_52_942_2_cropped.mp4';

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
tbytDat = tbytDat.data;
clearvars behFolder

%% get video file path and csv file path
filePath_vid = findFolderWithString(file_path, '_vid_cropped');
if isempty(filePath_vid)
    warning("The path for cropped videos doesn't exist!")
else
    filePath_vidFile = fullfile(filePath_vid, video_name);
end

% read csv file with coordinates
[~, video_name_base, ~] = fileparts(video_name);
filePath_csv = findCSVFileWithString(filePath_vid, video_name_base);
if ~isempty(filePath_csv)
    fprintf('Found the csv file for cooridnates!\n');
end

tableDlc = readDlcCsv(filePath_csv);

%% get behavioral events relative to frames
% get the number of trial
expression = '_(\d+)_cropped';  % Regular expression to match an integer between underscores

matches = regexp(filename, expression, 'tokens');
numTr = str2double(matches{1}{1});

faceCam_pt = evtInS.faceCam(evtInS.faceCam(:, 2)==numTr, 1);

vReader = VideoReader(filePath_vidFile);

if length(faceCam_pt) > vReader.NumFrames
    faceCam_pt = faceCam_pt(1:vReader.NumFrames);
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
% sprintf('%.2f', reFaceCam_pt(1))

%% generate the labeled video
vReader = VideoReader(filePath_vidFile);

labeledVideo = VideoWriter(fullfile(filePath_vid, strcat(video_name_base, '_', 'labeled')), 'MPEG-4');
labeledVideo.FrameRate = 200;
open(labeledVideo);

assert(vReader.NumFrames == size(tableDlc, 1));

dotColors = {[120, 95, 170], [81, 121, 235], [132, 171, 181], [149, 244, 207], ...
    [192, 249, 171], [231, 203, 128], [230, 131, 75], [206, 57, 32]}; % from the left clockwise
dotColors = cellfun(@(a) a./255, dotColors, 'un', 0);

% T


for ff = 1:vReader.NumFrames
    frame = readFrame(vReader);

    % draw the frame
    imshow(frame); hold on;

    dots = get_pupil_dots(tableDlc(ff, :), likelihoodcut);
    dots = cell2mat(dots);

    % draw dots
    for dt = 1:size(dots, 1)
        if ~isnan(dots(dt, 1))
            % Add dots to the frame with different colors and transparency
            scatter(dots(dt, 1), dots(dt, 2), 150, 'MarkerFaceColor', dotColors{dt}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5);
        end
    end
    hold off;

    % Capture the current figure as an image
    frameLabeled = getframe(gcf);
    frameLabeled = frameLabeled.cdata;

    % write the labeled frame to the output video
    writeVideo(labeledVideo, frameLabeled);
    fprintf('Frame #%d is completed.\n', ff);
end

% close the output video
close(labeledVideo);

%% get the pupil contour









end