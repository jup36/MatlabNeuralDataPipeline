function jsVideoPreprocess(filePath, nTr, varargin)

%filePath = 'C:\Users\parkj\Documents\DATA\VideoAPT\WR25_111618';
cd(filePath)

pVP = parse_input_jsVideo(filePath, nTr, varargin); 
% pVP = parse_input_jsVideo(filePath, 1, {'slowPlay',2}); 

% load the structure with video file info and Js kinematics (jsTime1k_KV)
if exist('jsTime1k_KV','var') ~= 1
    pathJsTime1k_KV = dir('**/*_Kinematics_VideoFiles.mat');
    S=load(fullfile(pathJsTime1k_KV.folder,pathJsTime1k_KV.name),'jsTime1k_KV');
    S=S.('jsTime1k_KV');    
else
    S=jsTime1k_KV; 
end

close all; 

% quick check on video data availability
if unique(isnan(S(nTr).fVideo)) || unique(isnan(S(nTr).sVideo)) || unique(isnan(S(nTr).vFrameTime)) || ~ismember(nTr,1:size(S,1))
    error('video files/data all missing!');
else
    fvid = VideoReader(S(nTr).fVideo);
    svid = VideoReader(S(nTr).sVideo);
end

% video file save dir
tbytVideoPath = fullfile(filePath,'tbytVideos');
if ~isfolder(fullfile(filePath,'tbytVideos'))
    mkdir(fullfile(filePath,'tbytVideos'))
end

% video file save name
vFileDateExp = '(20\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)'; % the expression of the date in the video file names
vFileDate = regexp(fvid.name,vFileDateExp,'match'); 
tbytVideoName = strcat(sprintf('Trial#%d',nTr),'_',S(nTr).trialType,'_',vFileDate{:});    

%videoPlayer = vision.VideoPlayer;
% new video
outputVideo = VideoWriter(fullfile(tbytVideoPath,tbytVideoName));
outputVideo.FrameRate = 250/pVP.Results.slowPlay; %fvid.FrameRate;
open(outputVideo);

% get the frame-by-frame time info
if strcmp(S(nTr).trialType,'sp') % for a successful trial, align time to pullStart
    frameTime = S(nTr).vFrameTime-(S(nTr).trJsReady+S(nTr).movKins.pullStart); % frameTime aligned to pullStart
    [~,pullStartI] = find(abs(frameTime)==min(abs(frameTime)));
    frameTimeC = num2cell(frameTime);
    frameTimeC = cellfun(@num2str,frameTimeC,'UniformOutput',false);
    frameTimeC = cellfun(@(c)strcat(c,'ms'),frameTimeC,'UniformOutput',false);
    [frameTimeC{pullStartI}] = deal('ReachStart');
    
elseif strcmp(S(nTr).trialType,'ps') % for a push trial, align time to pushStart
    frameTime = S(nTr).vFrameTime-(S(nTr).trJsReady+S(nTr).movKins.pushStart); % frameTime aligned to pushStart
elseif strcmp(S(nTr).trialType,'pmpp')  
    frameTime = S(nTr).vFrameTime-(S(nTr).trJsReady+S(nTr).movKins.pull.startI); % frameTime aligned to pullStart 
else  
    frameTime = S(nTr).vFrameTime-(S(nTr).trJsReady); % frameTime aligned to pullStart 
end
frameTimeC = num2cell(frameTime); 
frameTimeC = cellfun(@num2str,frameTimeC,'UniformOutput',false);
frameTimeC = cellfun(@(c)strcat(c,'ms'),frameTimeC,'UniformOutput',false);


txtPos = round([S(nTr).fVideoInfo.width/2 S(nTr).fVideoInfo.height/4]); % text position on frames

fr = 0; 
while hasFrame(fvid) && hasFrame(svid)
    fr = fr+1;
    img1 = readFrame(fvid);
    img2 = readFrame(svid);

    imgt = horzcat(img1, img2);
    
    RGB = insertText(imgt,txtPos,text_str,'FontSize',18,'BoxColor',...
    box_color,'BoxOpacity',0.4,'TextColor','white');
    
    % play video
    %step(videoPlayer, imgt);

    % record new video
    writeVideo(outputVideo, imgt);
end

% release
close(outputVideo);

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_jsVideo(filePath, nTr, vargs)
        % parse input, and extract name-value pairs
        default_frameRate = 250; % the default frame rate of videos
        default_slowPlay = 1; % fold to be slowed down for playback
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addRequired(p,'nTr')
        addParameter(p,'frameRate',default_frameRate)
        addParameter(p,'slowPlay',default_slowPlay)
        
        parse(p,filePath,nTr,vargs{:})    
    end
end