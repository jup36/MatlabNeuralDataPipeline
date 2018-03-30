function p = parse_input_beh( filePath, varargin )
% parse input, and extract name-value pairs

    p = inputParser; % create parser object
   
    default_numbNeuralProbe = 0;  % how many NIboard probes were used (e.g. zero if no NI neural probe was used)
    default_numbChEachProbe = 64; % how many channels on each NIboard probe
    default_XposCh = 33; % channel # for X position (default channel numbers for 64 channel recording)
    default_YposCh = 37; % channel # for Y position
    default_soleCh = 3;  % channel # for solenoid (water reward delivery)
    default_lickCh = 5;  % channel # for lick port
    default_laserCh = 7; % channel # for laser (laser TTL)
    default_numbTagLasers = 30; % the number of tagging trials given at the end of the experiment
    default_artifactRmv = true; % if true, removes the solenoid artifact from Xpos and Ypos channels by template subtraction
    
    addRequired(p,'filePath'); 
    addParameter(p,'numbNeuralProbe',default_numbNeuralProbe)
    addParameter(p,'numbChEachProbe',default_numbChEachProbe)
    addParameter(p,'XposCh',default_XposCh)
    addParameter(p,'YposCh',default_YposCh)
    addParameter(p,'soleCh',default_soleCh)
    addParameter(p,'lickCh',default_lickCh)
    addParameter(p,'laserCh',default_laserCh)
    addParameter(p,'numbTagLasers',default_numbTagLasers)
    addParameter(p,'artifactRmv',default_artifactRmv)
    
    parse(p,filePath,varargin{:})
    
end

