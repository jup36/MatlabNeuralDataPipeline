function TNC_GridCellReport(spikeMap,positions,rows,cols,thisRow);
% FUNCTION DETAILS:% simple function meant to generate a report of assorted properties from the spikeMap structure that quantify the properties of the resulting response map.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% 


figure(10); hold on;
i=((thisRow-1)*cols)+1;

    subplot(rows,cols,i:i+1);
    plot(positions.x,positions.y,'k',spikeMap.spikes.x,spikeMap.spikes.y,'r.','MarkerSize',16);
    axis([-30 200 -30 200])
    axis off
    if thisRow==1
        title('Spike Locations');
    end
    
i=i+2
    subplot(rows,cols,i);
    plot(spikeMap.spikes.isiH./length(spikeMap.spikes.isi),'k','LineWidth',1.5)
    axis tight
    axis([5 150 0 0.1]);
    if thisRow==1
        title('Interspike Intervals');
    end
    
    if thisRow==rows
        xlabel('ISI (ms)');
    end
    
i=i+1
    subplot(rows,cols,i:i+1);
    surf(spikeMap.smth.freq)
    shading flat
    axis([0 175 0 175 -1 15])
    view(0,90)
    axis off
    colormap(jet)
    if thisRow==1
        title('Frequency Map');
    end
    
i=i+2
    subplot(rows,cols,i);
    plot(spikeMap.smth.peakWidthsFreqHist./length(spikeMap.smth.peakWidthsFreqHist),'k','LineWidth',1.5)
    hold on
    plot(spikeMap.smth.fundSpaceDistFreqHist./length(spikeMap.smth.fundSpaceDistFreqHist),'r','LineWidth',1.5)
    axis([0 120 0 0.2])
    if thisRow==1
        title('Peak Properties');
    end
    
    if thisRow==rows
        xlabel('Distance (cm)');
    end
    
i=i+1
    subplot(rows,cols,i:i+1);
    surf(spikeMap.smth.auto)
    shading flat
    axis([0 350 0 350 -1 1])
    view(0,90)
    axis off
    colormap(jet)
    if thisRow==1
        title('Autocorrelation');
    end
    
i=i+2
    subplot(rows,cols,i);
    plot(spikeMap.smth.circSampAuto(2,:),spikeMap.smth.circSampAuto(1,:),'k-','LineWidth',1.5);
    axis([-pi pi -1 1]);
    if thisRow==1
        title('Circular correlation scores');
    end
    
    if thisRow==rows
        xlabel('Angle (rad)');
    end



