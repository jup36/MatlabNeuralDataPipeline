function p = parse_input_PP( filePath, vargs )
% parse input, and extract name-value pairs
default_reReadBin = false; % by default do not re-read the raw bin file, if done already
default_numbNeuralProbe = 0;  % specify how many NIboard probes were used (e.g. zero if no NI neural probe was used)
default_numbChEachProbe = 64; % specify how many channels are on the probe
default_cmosExpCh   = 1;    % ai0: CMOS exposure pulses
default_cmosTrigCh  = 2;    % ai1: CMOS trigger
default_speakerCh   = 4;    % ai3: reward delivery
default_lickCh      = 5;    % ai4: lick detector1
default_lick2Ch     = 7;    % ai6: lick detector2
default_bodyCamCh   = 6;    % ai5: bodyCam
default_faceCamCh   = 8;    % ai7: faceCam
default_photoDiodeCh = 15;  % ai14: photo diode
default_digitCh     = 17;   % the digital channel is always at the tailend

p = inputParser; % create parser object
addRequired(p,'filePath');
addParameter(p,'reReadBin',default_reReadBin);
addParameter(p,'numbNeuralProbe',default_numbNeuralProbe);
addParameter(p,'numbChEachProbe',default_numbChEachProbe);
addParameter(p,'cmosExpCh',default_cmosExpCh);
addParameter(p,'cmosTrigCh',default_cmosTrigCh);
addParameter(p,'speakerCh',default_speakerCh);
addParameter(p,'lickCh',default_lickCh);
addParameter(p,'lick2Ch',default_lick2Ch);
addParameter(p,'bodyCamCh',default_bodyCamCh);
addParameter(p,'faceCamCh',default_faceCamCh);
addParameter(p,'photoDiodeCh',default_photoDiodeCh);
addParameter(p,'digitCh',default_digitCh);

parse(p,filePath,vargs{:})
end