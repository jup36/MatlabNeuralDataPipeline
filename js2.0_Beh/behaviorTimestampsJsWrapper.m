function behaviorTimestampsJsWrapper(filePath, varargin)
%This is a wrapper function to run behaviorTimestampsJs.

defaultPath = 'Z:\parkj\NeuralData\js2.0';

if iscell(filePath)
    if isempty(filePath)
       filePath = {uigetdir(defaultPath)};         
    end
else
    error('Input a cell array with a filePath(s) or an empty cell array!'); 
end

for f = 1:length(filePath) 
    if isfolder(fullfile(filePath{f}))
        currentPath = fullfile(filePath{f}); 
    else
        error('Input a proper file path'); 
    end
    p = parse_input_Js(currentPath, varargin); 
    % p = parse_input_Js(currentPath, {'trialTimeout',10000}); 
    behaviorTimestampsJs(p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

    function p = parse_input_Js( filePath, vargs )
        % parse input, and extract name-value pairs
        default_numbNeuralProbe = 0;  % specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
        default_numbChEachProbe = 64; % specify how many channels are on the probe
        default_trStartCh = 33; % ch# for trial start
        default_camTrigCh = 34; % ch# for camera trigger
        default_rewardCh  = 35; % ch# for reward delivery
        default_trEndCh   = 36; % ch# for trial end
        default_encodeACh = 37; % ch# for stepper encoder A
        default_encodeBCh = 39; % ch# for stepper encoder B
        default_lickCh    = 1;  % ch# for lick detect (unattenuated channel)
        default_sgfiltFramelen = 101; % frame length for the sgolayfilt
        default_trialTimeout = 10000; % trial timeout duration
        default_pushThreshold = 50; % pushThreshold
        
        
        default_meanMass = [10 20 30 40 50 60 70 80 90 100; 2.95 3.52 4.10 4.67 5.25 5.82 6.40 6.97 7.55 8.12]; % torque(%) mass(g) mapping 
        
        p = inputParser; % create parser object
        addRequired(p,'filePath');
        addParameter(p,'numbNeuralProbe',default_numbNeuralProbe)
        addParameter(p,'numbChEachProbe',default_numbChEachProbe)
        addParameter(p,'trStartCh',default_trStartCh)
        addParameter(p,'camTrigCh',default_camTrigCh)
        addParameter(p,'rewardCh',default_rewardCh)
        addParameter(p,'trEndCh',default_trEndCh)
        addParameter(p,'encodeACh',default_encodeACh)
        addParameter(p,'encodeBCh',default_encodeBCh)
        addParameter(p,'lickCh',default_lickCh)
        addParameter(p,'sgfiltFramelen',default_sgfiltFramelen)
        addParameter(p,'trialTimeout',default_trialTimeout)
        addParameter(p,'pushThreshold', default_pushThreshold)
        addParameter(p,'meanMass', default_meanMass)
        
        parse(p,filePath,vargs{:})
    end


end
