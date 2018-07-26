function [] = TNC_DisplayExampleWaveforms(nevDS,numSpks)

colors = ['k' 'm' 'r' 'b' 'g' 'y' 'c'];

eArray = double(nevDS.Data.Spikes.Electrode);
uArray = double(nevDS.Data.Spikes.Unit);
wArray = double(nevDS.Data.Spikes.Waveform);
tArray = double(nevDS.Data.Spikes.Timestamps);

nnInds   = find(uArray>0);
allChan  = unique(eArray);
allTrode = allChan(find(allChan<129))
numSamps = size(wArray,2);   

dim      = ceil(sqrt(numel(allTrode)));
SpkPnts  = size(wArray,2)

winWidth = 30000.*10;
yScale   = 1000;
xScale   = 50;
startTime= tArray(1);

offsetX=0;
offsetY=0;

figure(1); clf;

disp(['Total recording duration: ' num2str((nevDS.Data.Spikes.Timestamps(numel(uArray))-nevDS.Data.Spikes.Timestamps(1))./30000./60) ' min.']);

for index = 1:numSpks

    if index==1
        for j=1:max(allTrode)
            plot([0 winWidth],[yScale.*j yScale.*j],'Color',[0.5 0.5 0.5]);   
            text(10,(yScale.*j)+(yScale./5),num2str(j));
            hold on;
            axis tight;
            axis off;            
        end
    end
    
    i = nnInds(index);

    title((nevDS.Data.Spikes.Timestamps(i)-nevDS.Data.Spikes.Timestamps(1))./30000);
    
    if eArray(i)<=64
        if uArray(i)<=numel(colors)
            offsetY = yScale.*eArray(i).*ones(1,SpkPnts);
            offsetX = tArray(i)-startTime;
            plot((1:SpkPnts)+offsetX,wArray(i,:)+offsetY,colors(uArray(i)));
        end
    end
    
    if rem(offsetX,10)==0
        drawnow;
        axis tight;
    end
    
    if offsetX>winWidth
        drawnow;
        startTime = tArray(i);
        pause();
        clf;
        for j=1:max(allTrode)
            plot([0 winWidth],[yScale.*j yScale.*j],'Color',[0.5 0.5 0.5]);            
            text(10,(yScale.*j)+(yScale./5),num2str(j));
            hold on;
            axis tight;
            axis off;
        end
    end
end