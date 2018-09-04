function p = parse_input_corrNeuralTrajMovKinematics( filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
%parse input, and extract name-value pairs for the main function 'corrNeuralTrajMovKinematics.m'

        default_trialFolds = 4;  % # of p.Results.trialFolds to classify the total trials by the rank of a behavioral variable
        default_lineWidth  = 2;  % default lineWidth to be used to plot the neural trajectories
        default_markerSize = 10; % default markerSize to be used to plot the neural trajectories
        default_saveNTjcFigs = false; % By default, there's no need to save the 3-d nTrj figures
        default_alpha = 0.05; % default alpha value to just plot the significant correlations between neural population trajectories and the behavioral kinematics
        
        p = inputParser; % create parser object
        
        addRequired(p,'filePath'); % file directory
        addRequired(p,'fileNameNeuralTrj'); % fileName for neural population trajectories
        addRequired(p,'fileNameBeh'); % fileName for movement kinematic variables
        addRequired(p,'saveNameTag'); % saveName used to save the outcomes
        addRequired(p,'eventMarkersRelative'); % timing for specific events to be marked on each neural trajectory
        addRequired(p,'reachTimeWin'); % timeWin in which the movement kinematics are concerned
        addRequired(p,'lickTimeWin');  % timeWin in which lick counts are concerned
        
        addParameter(p,'trialFolds', default_trialFolds) % the # of p.Results.trialFolds to divide the whole trials into
        addParameter(p,'lineWidth', default_lineWidth)   % the lineWidth to be used to plot nTrjs
        addParameter(p,'markerSize', default_markerSize) % the markerSzie to be used to plot nTrjs
        addParameter(p,'saveNtrjFigs', default_saveNTjcFigs) % By default, there's no need to save the 3-d nTrj figures
        addParameter(p,'alpha', default_alpha) % By default, the alpha level of 0.01 is used as the criterion for significant correlations
        
        parse(p, filePath, fileNameNeuralTrj, fileNameBeh, saveNameTag, eventMarkersRelative, reachTimeWin, lickTimeWin, vargs{:})

end