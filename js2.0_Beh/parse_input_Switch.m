function p = parse_input_Switch( filePath, vargs )
% parse input, and extract name-value pairs
default_reReadBin = false; % by default do not re-read the raw bin file, if done already
default_numbNeuralProbe = 0;  % specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
default_numbChEachProbe = 64; % specify how many channels are on the probe
default_camTrigCh = 1; % ch# for camera trigger
default_trStartCh = 2; % ch# for trial start (joystick ready)
default_encodeACh = 3; % ch# for stepper encoder A
default_encodeBCh = 4; % ch# for stepper encoder B
default_rewardCh  = 5; % ch# for reward delivery
default_trEndCh   = 6; % ch# for trial end
default_lickCh    = 7; % ch# for lick detect (unattenuated channel)
default_position1Ch = 8;
default_position2Ch = 9;
default_laserCh = 12;
default_pacerCh = 13;
default_sgfiltFramelen = 101; % frame length for the sgolayfilt
default_trialTimeout = 10000; % trial timeout duration
default_pushThreshold = 50; % pushThreshold
default_fps = 500; % default camera frame rate
default_laserUsed  = true; % laser used or not in the current experiment
default_numbTagLaser = 0; % the # of tagging lasers
default_tagLaserUsed = false; % laser tagging trials run in the current experiment
default_rewardDelay = 1000; % the reward TTL pulse is delivered right away not reflecting the delay, so add this to correct for it

default_meanMass = [10 20 30 40 50 60 70 80 90 100; 2.95 3.52 4.10 4.67 5.25 5.82 6.40 6.97 7.55 8.12]; % torque(%) mass(g) mapping

p = inputParser; % create parser object
addRequired(p,'filePath');
addParameter(p,'reReadBin',default_reReadBin);
addParameter(p,'numbNeuralProbe',default_numbNeuralProbe);
addParameter(p,'numbChEachProbe',default_numbChEachProbe);
addParameter(p,'camTrigCh',default_camTrigCh);
addParameter(p,'trStartCh',default_trStartCh);
addParameter(p,'encodeACh',default_encodeACh);
addParameter(p,'encodeBCh',default_encodeBCh);
addParameter(p,'rewardCh',default_rewardCh);
addParameter(p,'trEndCh',default_trEndCh);
addParameter(p,'lickCh',default_lickCh);
addParameter(p,'position1Ch', default_position1Ch);
addParameter(p,'position2Ch', default_position2Ch);
addParameter(p,'laserCh',default_laserCh);
addParameter(p,'pacerCh',default_pacerCh);
addParameter(p,'laserUsed',default_laserUsed);
addParameter(p,'tagLaserUsed',default_tagLaserUsed);
addParameter(p,'numbTagLaser',default_numbTagLaser);
addParameter(p,'sgfiltFramelen',default_sgfiltFramelen);
addParameter(p,'trialTimeout',default_trialTimeout);
addParameter(p,'pushThreshold',default_pushThreshold);
addParameter(p,'meanMass',default_meanMass);
addParameter(p,'fps',default_fps);
addParameter(p,'rewardDelay',default_rewardDelay)

parse(p,filePath,vargs{:})
end
