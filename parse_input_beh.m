function p = parse_input_beh( filePath, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function 
        % parse input, and extract name-value pairs
                
        default_numbNeuralProbe = 0;  % Specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
        default_numbChEachProbe = 64; % Specify the number of sites on the NIboard probe
        default_XposCh = 33; % channel # for X position (default channel numbers for 64 channel recording)
        default_YposCh = 37; % channel # for Y position
        default_pseudoLaserCh = 39; % channel # for pseudoLaser Pulses (without actual laser delivery)
        default_soleCh = 3;  % channel # for solenoid (water reward delivery)
        default_lickCh = 5;  % channel # for lick port
        default_laserCh = 7; % channel # for laser (laser TTL)
        default_camTrigCh = 40; % channel # for camera trigger 
        default_numbTagLasers = 60; % the number of tagging trials given at the end of the experiment
        default_artifactRmv = true; % if true, removes the solenoid artifact from Xpos and Ypos channels by template subtraction
        default_reachBeforeLastReward = true; % logical to detect reaches before the last reward delivery
        default_laserUsed = true; % logical to indicate whether laser stimulation was used during the session or not
        default_highPassFilterFC = 400; % 300Hz e.g., a highPass filter can be applied to denoise a channel
        default_filterLickCh = false; % logical to apply a highPass filter to the lick ch
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'numbNeuralProbe',default_numbNeuralProbe)
        addParameter(p,'numbChEachProbe',default_numbChEachProbe)
        addParameter(p,'XposCh',default_XposCh)
        addParameter(p,'YposCh',default_YposCh)
        addParameter(p,'pseudoLaserCh', default_pseudoLaserCh)
        addParameter(p,'soleCh',default_soleCh)
        addParameter(p,'lickCh',default_lickCh)
        addParameter(p,'laserCh',default_laserCh)
        addParameter(p,'camTrigCh',default_camTrigCh)
        addParameter(p,'numbTagLasers',default_numbTagLasers)
        addParameter(p,'artifactRmv',default_artifactRmv)
        addParameter(p,'reachBeforeLastReward',default_reachBeforeLastReward)
        addParameter(p,'laserUsed', default_laserUsed)
        addParameter(p,'highPassFilterFC', default_highPassFilterFC)
        addParameter(p,'filterLickCh', default_filterLickCh)
        
        parse(p,filePath,vargs{:})
        
    end