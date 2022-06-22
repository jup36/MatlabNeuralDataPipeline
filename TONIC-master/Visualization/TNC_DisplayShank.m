function [] = TNC_DisplayShank(fileNameStr,arrayType,shankNumber,color,timeRange,downsample,filter,yscaling,figNum)
% FUNCTION DETAILS: Plots data from all electrodes on a given shank to the current axes
%     fileNameStr >> name of the ns5 file from which data is loaded
%     arrayType >> electrode array type: 'NN_w64', 'NN_b64', 'NN_w32', (see TNC_RemapElecPos for details)
%     shankNumber >> Which shank to plot?
%     color >> color for display
%     pntRange >> number of points to display
%     downsample >> decimation factor for the display
%     filter >> high pass filter the data == 1, broadband == 0
%     yscaling >> amount to scale in y between channels (if <10 assumed to be number of s.d.s if >10 assumed to be in uV)

[row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);
if shankNumber(2)==0
    eArray = electrodeMatrix(:,shankNumber(1));
    disp(['Plotting data for shank ' num2str(shankNumber(1)) '...']);
else
    eArray = electrodeMatrix(shankNumber(1),:);
    disp(['Plotting data for electrode row ' num2str(shankNumber(1)) '...']);
end

figure(figNum); hold off;

for i=2:numel(eArray)

    % LOAD DATA INTO THE Ns5DATA STRUCTURE    
    disp(['Loading data from channel: ' num2str(eArray(i))]);
    disp(['...over the time range: ' 't: ' num2str(timeRange(1)) ':' num2str(timeRange(2))]);
    Ns5DATA = openNSx('report','read',fileNameStr,['e: ' num2str(eArray(i))],['t: ' num2str(timeRange(1)) ':' num2str(timeRange(2))],'min');
    rawData = sgolayfilt(Ns5DATA.Data(1,:),11,21);
    forSpec.data(i-1,:) = rawData;
    
    if filter==1
        % design the bandpass filter to use for field potentials
        dLo = fdesign.bandpass('n,fc1,fc2', 2000, 22, 32, Ns5DATA.MetaTags.SamplingFreq);
        HdLo = design(dLo);
        betaBandData.lowCutOff = 22
        betaBandData.hiCutOff = 32
        betaBandData.values = filtfilt(HdLo.Numerator,1,rawData); %zero-phase filtering
        dispData(i,:) = decimate(betaBandData.values,downsample);
    else
        dispData(i,:) = decimate(rawData,downsample);
    end
    
    if yscaling<10
        yscaling = yscaling;
    end
    
    figure(figNum); 
    subplot(numel(eArray)+2,1,1:numel(eArray));
    plot(1:size(dispData(1,:),2),(dispData(i,:)./std(dispData(i,:)))+(yscaling*(i-2)),'Color',color); hold on;
    title(['e: ' num2str(eArray(i))]);
    axis tight; drawnow;
    
end

figure(figNum+10)
imagesc(corr(dispData'),[0 1]);

figure(figNum+20)
wname  = 'cgau2';
scales = 1:53;
ntw = 21;
wcoher(dispData(2,:),dispData(7,:),scales,wname,'ntw',ntw,'nsw',1,'plot','wcs');

figure(figNum+30)
paramsC.tapers      = [7 13];
paramsC.pad         = 0;
paramsC.Fs          = 500; % reflects the fact that the data was decimated down to 1 kHz sampling frequency
paramsC.fpass       = [0 100];
paramsC.err         = [1 0.05];
% movingwin           = [0.4 0.01];

for i=1:size(forSpec.data,1)
    [S,f,Serr]=mtspectrumc(decimate(forSpec.data(i,:),60),paramsC);
    forSpec.spec(i,:) = S;
end
semilogy(f,mean(forSpec.spec)); axis tight;

% 
% % LOAD the joystick data
% Ns5DATA = openNSx('report','read',fileNameStr,'e: 138',['t: ' num2str(timeRange(1)) ':' num2str(timeRange(2))],'min');
% rawData = sgolayfilt(Ns5DATA.Data(1,:),11,21);
% levData(1,:) = decimate(rawData,downsample);
% 
% subplot(numel(eArray)+1,1,i+1);
% plot(1:size(levData,2),levData(1,:),'Color',color); hold on;
% title(['e: 137']);
% axis tight; drawnow;
