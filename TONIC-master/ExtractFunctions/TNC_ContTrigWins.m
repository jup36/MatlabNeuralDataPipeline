function [sink] = TNC_ContTrigWins(fileNameStr,times,window,chanList,ContData)

% If you want to just use the center 6 shanks use 
%     chanList = electrodeMap(9:56);

numTimes            = numel(times);
numTrodes           = numel(chanList);
paramsC.tapers      = [5 9];
paramsC.pad         = 0;
paramsC.Fs          = 1000; % reflects the fact that the data was decimated down to 1 kHz sampling frequency
paramsC.fpass       = [2 115];
paramsC.err         = 0;
movingwin           = [0.2 0.01];

if numTimes==1
    dispFlag=1;
else
    dispFlag=0;
end

SamplingFreq = 30000;
sBT = fdesign.bandpass('n,fc1,fc2', 1000, 250, 7500, SamplingFreq);
HsBT = design(sBT);
dBT = fdesign.bandpass('n,fc1,fc2', 300, 18, 32, SamplingFreq./30);
HdBT = design(dBT);
dGL = fdesign.bandpass('n,fc1,fc2', 300, 40, 55, SamplingFreq./30);
HdGL = design(dGL); 
dGH = fdesign.bandpass('n,fc1,fc2', 300, 65, 115, SamplingFreq./30);
HdGH = design(dGH);  

        
for i=1:numTimes
    
    currTime = times(i);
    begT = round(currTime-window(1));
    endT = round(currTime+window(2));
    
    % Load the data over the -window(1) to window(2) time frame
    timeRange = ['t: ' sprintf('%8.4f',begT./1000) ':' sprintf('%8.4f',endT./1000)];
    disp(['Reach ' num2str(i) ' of ' num2str(numTimes) ' | at ' timeRange]);
    Ns5DATA = openNSx('read',fileNameStr, timeRange, 'sec'); 
    if i==1
        minSamps = ((window(1)+window(2)).*round(Ns5DATA.MetaTags.SamplingFreq./1000))-1
    end

    % retrieve licks and voltage traces just to double check
    sink.behavior.reach(i).vel  = ContData.behavior.sLeverV(begT:endT);
    sink.behavior.lick(i).raw   = ContData.behavior.rawLick(begT:endT);

    % denoise by subtracting the mode of all channels?
    avgSignal = mean(Ns5DATA.Data(:,1:minSamps),1);
    
    for k=1:numTrodes

        j=chanList(k);
        
        % Filter tight in the beta range
        spkDataTMP = filtfilt(HsBT.Numerator,1,Ns5DATA.Data(j,1:minSamps)-avgSignal); % zero-phase filtering
        spkData = sgolayfilt(spkDataTMP,11,21);
        
        % Filter to show MUA and lowband
        lfpData = decimate(Ns5DATA.Data(j,1:minSamps),30);

        % Decimate lowband and calculate spectrogram
        [S,t,f] = mtspecgramc(lfpData,movingwin,paramsC);

        % Filter tight in the beta range
        beta = filtfilt(HdBT.Numerator,1,lfpData); % zero-phase filtering
        gammaHi = filtfilt(HdGH.Numerator,1,lfpData); % zero-phase filtering  
        gammaLo = filtfilt(HdGL.Numerator,1,lfpData); % zero-phase filtering

        sink.phys.trode(j).beta(i,:)    = beta;
        sink.phys.trode(j).glo(i,:)     = gammaLo;
        sink.phys.trode(j).ghi(i,:)     = gammaHi;
        sink.phys.trode(j).spkData(i,:) = spkData;
        sink.phys.trode(j).spec(i).S    = S;
        sink.phys.trode(j).spec(i).t    = t;
        sink.phys.trode(j).spec(i).f    = f;
        
        % Get estimates of instantaneous power by using a median filter of squared amplitudes
        % uses a window of 80 points (3-4 cycles of beta frequency), or 0.25 sec at 500 Hz sampling rate
        sink.phys.trode(j).betaP(i,:)   = medfilt1(beta.^2,80);
        sink.phys.trode(j).gloP(i,:)    = medfilt1(gammaLo.^2,40);
        sink.phys.trode(j).ghiP(i,:)    = medfilt1(gammaHi.^2,25);        
    
        % optional plotting
        if dispFlag==1
            figure(2)
            subplot(611);
                plot((-window(1)).*1000:(+window(2)).*1000,sink.behavior.reach(i).vel); axis tight;
                title(['Electrode: ' num2str(j)]);
                ylabel('Velocity (a.u.)');
                xlabel('Time (ms)');
            subplot(612);
                plot((-window(1)).*1000:(+window(2)).*1000,sink.behavior.lick(i).raw); axis tight;
                ylabel('Lickometer (uV)');
                xlabel('Time (ms)');
            subplot(613);
                plot(sink.phys.trode(j).spkData(i,:)); axis([0 6e4 -750 250]);
                ylabel('Voltage (uV)');
                xlabel('Samples (30 kHz resolution)');
            subplot(614);
                imagesc((sink.phys.trode(j).spec(i).t-1).*1000,sink.phys.trode(j).spec(i).f,log10(sink.phys.trode(j).spec(i).S')); axis tight;
                ylabel('Frequency (Hz)');
                xlabel('Time (ms)');
            subplot(615);
                plot((-window(1)).*1000:(+window(2)).*1000-1,sink.phys.trode(j).betaP(i,:),'r'); axis([-1e3 2e3 0 13000]);
                ylabel('Beta Power (mV^2)');
                xlabel('Time (ms)');
            subplot(616);
                plot((-window(1)).*1000:(+window(2)).*1000-1,sink.phys.trode(j).ghiP(i,:),'k'); axis([-1e3 2e3 0 3.5e3]);
                ylabel('Gamma Power (mV^2)');
                xlabel('Time (ms)');
            drawnow; pause();
        end

    end
                
end
