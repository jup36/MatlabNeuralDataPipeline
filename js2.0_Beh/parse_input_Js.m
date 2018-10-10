function p = parse_input_Js( filePath, vargs ) % note that a nested function must use vargs not varargin when varargin was used for the main function
% parse input, and extract name-value pairs
default_numbNeuralProbe = 0;  % specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
default_numbChEachProbe = 64; % specify how many channels are on the probe
default_trStartCh = 33; % ch# for trial start
default_rewardCh  = 35; % ch# for reward delivery
default_trEndCh   = 36; % ch# for trial end
default_dirCh     = 37; % ch# for stepper dir
default_stepCh    = 39; % ch# for stepper steps
default_lickCh    = 1;  % ch# for lick detect (unattenuated channel)

p = inputParser; % create parser object
addRequired(p,'filePath');
addParameter(p,'numbNeuralProbe',default_numbNeuralProbe)
addParameter(p,'numbChEachProbe',default_numbChEachProbe)
addParameter(p,'trStartCh',default_trStartCh)
addParameter(p,'rewardCh',default_rewardCh)
addParameter(p,'trEndCh',default_trEndCh)
addParameter(p,'dirCh',default_dirCh)
addParameter(p,'stepCh',default_stepCh)
addParameter(p,'lickCh',default_lickCh)

parse(p,filePath,vargs{:})
end