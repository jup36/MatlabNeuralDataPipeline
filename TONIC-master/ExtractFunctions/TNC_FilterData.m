function [lowBandData,hiBandData] = TNC_FilterData(data,sampleRate,decimateFlag,plotFlag,bandFlag)
% FUNCTION DETAILS: general utility to separate data into a pair of bandwidths for spike sorting and continuous analysis. LowBand is 2-0.1k; HiBand: 0.7k-7k
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

tstart = tic;

% design the bandpass filter to use for spike data (use fvtool on the filter to visualize)
    if sampleRate < 30000
        dHi = fdesign.bandpass('n,fc1,fc2', 1000, 400, (sampleRate./2)-500, sampleRate);
        disp(['High pass cut-off set at ' num2str((sampleRate./2)-500) ' Hz.']);
    else
        dHi = fdesign.bandpass('n,fc1,fc2', 1000, 400, 7000, sampleRate);
    end
    HdHi = design(dHi);
    hiBandData.lowCutOff = 400;
    hiBandData.hiCutOff = 7000;
    hiBandData.sampleRate = sampleRate;

    if bandFlag(1)==1
        % design the bandpass filter to use for field potentials
        dLo = fdesign.bandpass('n,fc1,fc2', 2000, 2, 110, sampleRate);
        HdLo = design(dLo);
        lowBandData.lowCutOff = 2;
        lowBandData.hiCutOff = 110;
        hiBandData.values = filtfilt(HdHi.Numerator,1,data); %zero-phase filtering    
    else
        hiBandData.values = [];
        hiBandData.sampleRate = sampleRate; 
    end

    if bandFlag(2)==1
        % run the filter on the data supplied at the call
        lowBandDataTMP = filtfilt(HdLo.Numerator,1,data); %zero-phase filtering

        % after filtering decimate the field potential data to ~1kHz
        if decimateFlag>0
            lowBandData.values = decimate(lowBandDataTMP,decimateFlag);
            lowBandData.sampleRate = sampleRate./decimateFlag;
        else
            lowBandData.values = lowBandDataTMP;
            lowBandData.sampleRate = sampleRate;
        end
    else
        lowBandData.values = [];
        lowBandData.sampleRate = [];
    end

    hiBandData.sampleRate = sampleRate;
    
% plot the resulting data
    if plotFlag
        figure(1); clf;
        plot(1:1:length(data),data(1,:),'k',1:1:length(hiBandData.values),hiBandData.values(1,:),'b',1:length(lowBandData.values),lowBandData.values(1,:),'r')
    end

    telapsed = toc(tstart);

% report the status and elapsed time
    disp(sprintf('     EVENT FILTER | Data size (seconds): %g | Elapsed time (seconds): %g',length(data)./sampleRate,telapsed));
    
    
% %% Consider using wavelet filters in the future:
% % Matlab code for wavelet filtering.
% % This function requires the Wavelet Toolbox.
% 
% function fdata = wavefilter(data, maxlevel)
% % fdata = wavefilter(data, maxlevel)
% % data	- an N x M array of continuously-recorded raw data
% %		where N is the number of channels, each containing M samples
% % maxlevel - the level of decomposition to perform on the data. This integer
% %		implicitly defines the cutoff frequency of the filter.
% % 		Specifically, cutoff frequency = samplingrate/(2^(maxlevel+1))
% 
% [numwires, numpoints] = size(data);
% fdata = zeros(numwires, numpoints);
% 
% % We will be using the Daubechies(4) wavelet.
% % Other available wavelets can be found by typing 'wavenames'
% % into the Matlab console.
% wname = 'db4'; 
% 
% for i=1:numwires % For each wire
%     % Decompose the data
%     [c,l] = wavedec(data(i,:), maxlevel, wname);
%     % Zero out the approximation coefficients
%     c = wthcoef('a', c, l);
%     % then reconstruct the signal, which now lacks low-frequency components
%     fdata(i,:) = waverec(c, l, wname);
% end