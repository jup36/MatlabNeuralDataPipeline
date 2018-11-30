    function p = parse_input_jsVideo(filePath, numbTrial, vargs)
        % parse input, and extract name-value pairs
        default_frameRate = 250; % the default frame rate of videos
        default_slowPlay = 5; % fold to be slowed down for playback
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addRequired(p,'numbTrial')
        addParameter(p,'frameRate',default_frameRate)
        addParameter(p,'slowPlay',default_slowPlay)
        
        parse(p,filePath,numbTrial,vargs{:})    
    end