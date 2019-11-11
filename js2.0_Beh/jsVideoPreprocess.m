function jsVideoPreprocess(filePath, nTr, varargin)

%filePath = 'Z:\parkj\NeuralData\js2.0\WR25\110718_LowHighShift';
cd(filePath)

pVP = parse_input_jsVideo(filePath, nTr, varargin);
% pVP = parse_input_jsVideo(filePath, 19, {'slowPlay',10, 'frTimeRange', [-1000 1000]});

% load the structure with video file info and Js kinematics (jsTime1k_KV)
if exist('jsTime1k_KV','var') ~= 1
    pathJsTime1k_KV = dir('**/*_Kinematics_VideoFiles.mat');
    S=load(fullfile(pathJsTime1k_KV.folder,pathJsTime1k_KV.name),'jsTime1k_KV');
    S=S.('jsTime1k_KV');
else
    S=jsTime1k_KV;
end
close all;

% video file save dir
tbytVideoPath = fullfile(filePath,'tbytVideos');
if ~isfolder(fullfile(filePath,'tbytVideos'))
    mkdir(fullfile(filePath,'tbytVideos'))
end

for i = 1:length(nTr)
    % quick check on video data availability
    if unique(isnan(S(nTr(i)).fVideo)) || unique(isnan(S(nTr(i)).sVideo)) || unique(isnan(S(nTr(i)).vFrameTime)) || ~ismember(nTr(i),1:length(S))
        warning(strcat('some data are missing for trial#',num2str(nTr(i)),'-Video not produced!'))
        continue;
        %error('video files/data all missing!');
    else
        fvid = VideoReader(S(nTr(i)).fVideo);
        svid = VideoReader(S(nTr(i)).sVideo);
    end
        
    % video file save name
    vFileDateExp = '(20\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)'; % the expression of the date in the video file names
    vFileDate = regexp(fvid.name,vFileDateExp,'match');
    tbytVideoName = strcat(sprintf('Trial#%d',nTr(i)),'_',S(nTr(i)).trialType,'_',vFileDate{:});
    
    %videoPlayer = vision.VideoPlayer;
    % new video
    outputVideo = VideoWriter(fullfile(tbytVideoPath,tbytVideoName),'MPEG-4');
    outputVideo.FrameRate = 250/pVP.Results.slowPlay; %fvid.FrameRate;
    outputVideo.Quality = 100;
    open(outputVideo);
    
    if isempty(pVP.Results.frTimeRange)
    elseif length(pVP.Results.frTimeRange)==2
        frameTimeRange = linspace(pVP.Results.frTimeRange(1), pVP.Results.frTimeRange(2), pVP.Results.frTimeRange(2)-pVP.Results.frTimeRange(1)+1 );
    else
        error('To select frames, define the START & END of the frameTimeRange!')
    end
    
    % get the frame-by-frame time info
    if strcmp(S(nTr(i)).trialType,'sp') % for a successful trial, align time to pullStart
        frameTime = S(nTr(i)).vFrameTime-(S(nTr(i)).trJsReady+S(nTr(i)).movKins.pullStart); % frameTime aligned to pullStart
        pullTime = 0:S(nTr(i)).movKins.pullStop - S(nTr(i)).movKins.pullStart; % pull time bins aligned to pullStart
        pullTimeFrames = ismember(frameTime,pullTime); % logical indicating pullTime frames      
        % identify the frames to include
        frameTimeLogic = ones(1,length(frameTime)); % by default, include all frames
        if ~isempty(pVP.Results.frTimeRange)
            frameTimeLogic = ismember(frameTime,frameTimeRange);
        end
        [~,pullStartI] = find(abs(frameTime)==min(abs(frameTime))); % mark the frames closest to pullStart
        frameTimeC = num2cell(frameTime); % frameTime cell
        frameTimeC = cellfun(@num2str,frameTimeC,'UniformOutput',false);
        frameTimeC = cellfun(@(c)strcat(c,'ms'),frameTimeC,'UniformOutput',false);
        [frameTimeC{pullStartI}] = deal('ReachStart');
        
    elseif strcmp(S(nTr(i)).trialType,'ps') % for a push trial, align time to pushStart
        frameTime = S(nTr(i)).vFrameTime-(S(nTr(i)).trJsReady+S(nTr(i)).movKins.pushStart); % frameTime aligned to pushStart
        pushTime = 0:S(nTr(i)).movKins.pushStop - S(nTr(i)).movKins.pushStart; % push time bins aligned to pullStart
        pushTimeFrames = ismember(frameTime,pushTime); % logical indicating pullTime frames   
        % identify the frames to include
        frameTimeLogic = ones(1,length(frameTime)); % by default, include all frames
        if ~isempty(pVP.Results.frTimeRange)
            frameTimeLogic = ismember(frameTime,frameTimeRange);
        end
        [~,pushStartI] = find(abs(frameTime)==min(abs(frameTime))); % mark the frames closest to pushStart
        frameTimeC = num2cell(frameTime); % frameTime cell
        frameTimeC = cellfun(@num2str,frameTimeC,'UniformOutput',false);
        frameTimeC = cellfun(@(c)strcat(c,'ms'),frameTimeC,'UniformOutput',false);
        [frameTimeC{pushStartI}] = deal('pushStart');
        
    elseif strcmp(S(nTr(i)).trialType,'pmpp')
        frameTime = S(nTr(i)).vFrameTime-(S(nTr(i)).trJsReady+S(nTr(i)).movKins.pull.startI); % frameTime aligned to pullStart
        pullTime = 0:S(nTr(i)).movKins.pull.stopI - S(nTr(i)).movKins.pull.startI; % pull time bins aligned to pullStart
        pullTimeFrames = ismember(frameTime,pullTime); % logical indicating pullTime frames     
        [~,pullStartI] = find(abs(frameTime)==min(abs(frameTime))); % mark the frames closest to pullStart
        frameTimeC = num2cell(frameTime); % frameTime cell
        frameTimeC = cellfun(@num2str,frameTimeC,'UniformOutput',false);
        frameTimeC = cellfun(@(c)strcat(c,'ms'),frameTimeC,'UniformOutput',false);
        [frameTimeC{pullStartI}] = deal('ReachStart');
        % identify the frames to include
        frameTimeLogic = ones(1,length(frameTime)); % by default, include all frames
        if ~isempty(pVP.Results.frTimeRange)
            frameTimeLogic = ismember(frameTime,frameTimeRange);
        end
    else % time-out
        frameTime = S(nTr(i)).vFrameTime-(S(nTr(i)).trJsReady); % frameTime aligned to js ready time
        % identify the frames to include
        frameTimeLogic = ones(1,length(frameTime)); % by default, include all frames
        if ~isempty(pVP.Results.frTimeRange)
            frameTimeLogic = ismember(frameTime,frameTimeRange);
        end
        frameTimeC = num2cell(frameTime); % frameTime cell
        frameTimeC = cellfun(@num2str,frameTimeC,'UniformOutput',false);
        frameTimeC = cellfun(@(c)strcat(c,'ms'),frameTimeC,'UniformOutput',false);
    end
    
    % determine at which frame to start recording
    if unique([S(nTr(i)).vFronFileCalled,S(nTr(i)).vSideFileCalled])==2 % if the video file got called previously, align the frame to the beginning of the relative frames among all frames of the video file
        fvid.CurrentTime = (S(nTr(i)).origVideoFramesCnt-length(frameTimeC))+1;
        svid.CurrentTime = fvid.CurrentTime; 
    else % in most cases just start recording from the beginning
        fvid.CurrentTime = 1;
        svid.CurrentTime = 1; 
    end
    
    txtPos = round([S(nTr(i)).fVideoInfo.width*19/20 S(nTr(i)).fVideoInfo.height*9/10]); % text position on frames
    circlePos = round([S(nTr(i)).fVideoInfo.width*19/20 S(nTr(i)).fVideoInfo.height*1/15]); % circle position on frames to indicate frames of action 
    slowFactor = sprintf('%.2f',1/pVP.Results.slowPlay); % to display the playback speed up to two decimal points
    fr = 0;
    while hasFrame(fvid) && hasFrame(svid) && fr <= length(frameTimeC)
    %for fr = find(frameTimeLogic==1)
        fr = fr+1;
        
        img1 = readFrame(fvid);
        %img1 = mmread(fullfile(fvid.Path,fvid.name),fr);
        img2 = readFrame(svid);
        %img2 = mmread(fullfile(svid.Path,svid.name),fr);
        
        if frameTimeLogic(fr) % if the current frame's to be included
            imgt = horzcat(img1, img2);
            %imgt = horzcat(img1.frames.cdata, img2.frames.cdata); 
            
            if isempty(strfind(frameTimeC{fr},'ms'))
                text_color = 'Red';
            else
                text_color = 'Yellow';
            end
            
            imshow(imgt);
            text(1,txtPos(2),strcat(frameTimeC{fr},'_x',slowFactor),'FontSize',18,'Color',text_color,'Interpreter', 'none','LineStyle','none'); % insert text indicating time/playback speed info
            hold on; 
            % add a circle to indicate action frames
            if strcmp(S(nTr(i)).trialType,'sp')
                if pullTimeFrames(fr) 
                    vc = scatter(20,circlePos(2),300,[0 1 1],'filled'); 
                    alpha(vc,.5)
                end
            elseif strcmp(S(nTr(i)).trialType,'ps')
                if pushTimeFrames(fr)
                    vc = scatter(20,circlePos(2),300,[1 0 1],'filled'); 
                    alpha(vc,.5)
                end
            elseif strcmp(S(nTr(i)).trialType,'pmpp')
                if pullTimeFrames(fr)
                    vc = scatter(20,circlePos(2),300,[0 1 1],'filled'); 
                    alpha(vc,.5)
                end
            end
            hold off; 
            
            thisFrame = getframe(gca);
            % play video
            %step(videoPlayer, imgt);
            
            % record new video
            writeVideo(outputVideo, thisFrame);
        else
        end
    end
    % release
    close(outputVideo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
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
end