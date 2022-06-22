%% PARAMS
snrThresh = 3;
outputStructure.snrThresh = snrThresh;

%% DATA LOADING
% filenamestr = '/Volumes/huxley/m42-SomeData/m42-20110427-wide-plustim-008.ns5'
outputStructure.filenamestr = filenamestr;

%% EXTRACT MUA AND LFP DATA

for i=1:numel(chanList)
    
    % Load the data for a given channel
    disp('   ');
    disp('   ');

    disp(['Loading data from electrode ' num2str(i) '...']);
    electrode = ['e:' num2str(i)];
    dataSpks = openNSx('report','read',filenamestr, electrode);

    % Filter the data to select for spikes
    disp('Filtering data...');
    [lowBandData,hiBandData] = TNC_FilterData(dataSpks.Data(1,:),dataSpks.MetaTags.SamplingFreq,1,0);
    paramsC.tapers   = [3 5];
    paramsC.pad      = 0;
    paramsC.Fs       = lowBandData.sampleRate;
    paramsC.fpass    = [0 100];
    paramsC.err      = 0;
    movingwin        = [0.5 0.05];
    data = rmlinesmovingwinc(lowBandData.values,movingwin,10,paramsC,0.05,'n',60);
    lowBandData.values=data';

    % Find threshold crossings with a specific SNR
    disp('Finding threshold crossings...');
    [events] = TNC_EventDetect(hiBandData.values,hiBandData.sampleRate,snrThresh);

    % Store the delta function sampled at 1kHz and the lfp amplitude
    tmpTs = events.inds .* events.sampleRate;
    delta = zeros(1,round(size(dataSpks.Data(1,:),2).*(1000./events.sampleRate)));
    delta(1,round(events.inds.*(1000./events.sampleRate))) = 1;
    outputStructure.channel(i).delta = delta;

    % Apply chronux multitaper method to extract power spectrum in time
    disp('Calculating the spectrogram of the lowBand data...');
    movingwin       = [1 0.01];
    outputStructure.channel(i).lfp.params      = paramsC;
    outputStructure.channel(i).lfp.movingwin   = movingwin;
    [S,t,f] = mtspecgramc(lowBandData.values,movingwin,paramsC);

    figure(1);
    imagesc(t,f,S',[0 60]);

    % Store the amplitude window of specific freq bands as a function of time
    outputStructure.channel(i).lfp.t = t;
    outputStructure.channel(i).lfp.f = f;
    outputStructure.channel(i).lfp.S = S;

end

%% Load the event data to get trial boundaries and store
filenamestrE = [filenamestr(1,1:length(filenamestr)-3) 'nev']
dataEvents = openNEV(filenamestrE,'report','read','nomat');

nevRes = dataEvents.MetaTags.SampleRes;

rewardInds = find(dataEvents.Data.Spikes.Electrode==142); % reward
outputStructure.behavior.rewardInds = dataEvents.Data.Spikes.Timestamps(rewardInds).*(1000./nevRes);

threshInds = find(dataEvents.Data.Spikes.Electrode==143); % threshold crossing
outputStructure.behavior.threshInds = dataEvents.Data.Spikes.Timestamps(threshInds).*(1000./nevRes);

lickInds = find(dataEvents.Data.Spikes.Electrode==139); % licking
outputStructure.behavior.lickInds = dataEvents.Data.Spikes.Timestamps(lickInds).*(1000./nevRes);

%% Write the output structure
disp(['EVAL: save ' filenamestr(1,1:length(filenamestr)-4) '-MUAextract.mat outputStructure']);
eval(['save ' filenamestr(1,1:length(filenamestr)-4) '-MUAextract.mat outputStructure']);