function jsVideoPreprocessSimpleMerge(filePath,fVfilePath, sVfilePath, slowPlay)

%filePath = 'Z:\parkj\NeuralData\js2.0\WR25\110718_LowHighShift';
cd(filePath)

% video file save dir
tbytVideoPath = fullfile(filePath,'tbytVideos');
if ~isfolder(fullfile(filePath,'tbytVideos'))
    mkdir(fullfile(filePath,'tbytVideos'))
end


if exist(fullfile(fVfilePath),'file')==2 && exist(fullfile(sVfilePath),'file')==2
    fvid = VideoReader(fullfile(fVfilePath));
    svid = VideoReader(fullfile(sVfilePath));
    
    % video file save name
    vFileDateExp = '(20\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)'; % the expression of the date in the video file names
    vFileDate = regexp(fvid.name,vFileDateExp,'match');
    tbytVideoName = strcat('simpleMerge','_',vFileDate{:});
    
    % new video
    outputVideo = VideoWriter(fullfile(tbytVideoPath,tbytVideoName));
    outputVideo.FrameRate = 250/slowPlay; %fvid.FrameRate;
    outputVideo.Quality = 100;
    open(outputVideo);
    
    % determine at which frame to start recording
    fvid.CurrentTime = 1;
    svid.CurrentTime = 1;
    
    slowFactor = sprintf('%.2f',1/slowPlay); % to display the playback speed up to two decimal points
    fr = 0;
    timeElapsed = 0; 
    while hasFrame(fvid) && hasFrame(svid)
        fr = fr+1;
        timeElapsed = timeElapsed + 4; %1000/250  
        
        img1 = readFrame(fvid);
        %img1 = mmread(fullfile(fvid.Path,fvid.name),fr);
        img2 = readFrame(svid);
        %img2 = mmread(fullfile(svid.Path,svid.name),fr);
        
        imgt = horzcat(img1, img2);
        %imgt = horzcat(img1.frames.cdata, img2.frames.cdata);
        
        imshow(imgt);
        text(1,fvid.height*9/10,strcat(sprintf('%dms',timeElapsed),'_x',slowFactor),'FontSize',18,'Color',[255 211 0]./255,'Interpreter', 'none','LineStyle','none'); % insert text indicating time/playback speed info
        
        thisFrame = getframe(gca);
        % play video
        %step(videoPlayer, imgt);
        
        % record new video
        writeVideo(outputVideo, thisFrame);
    end
    % release
    close(outputVideo);
end