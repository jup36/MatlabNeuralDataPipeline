% plot a bunch of traces
numTraces = 64;
figure(1); clf;
figure(2); clf;

for i=1:64
    [lowBandData,hiBandData] = TNC_FilterData(NS5data.Data(i,:),0);
    lowbandDataTogether(i,:) = lowBandData.values(1,:);
    hibandDataTogether(i,:) = hiBandData.values(1,:);
%     figure(1);
%     subplot(64,1,i);
%     plot(1:60000,lowbandDataTogether(i,:));
%     axis tight
%     axis off
%     
%     figure(2);
%     subplot(64,1,i);
%     plot(1:60000,hibandDataTogether(i,:));
%     axis tight
%     axis off
    
    filenameh = sprintf('%s_hi.h5','data');
    name = sprintf('/RecordA%g',i);
    if(i==1)
        hdf5write(filenameh,name,hibandDataTogether(i,:));
    else
        hdf5write(filenameh,name,hibandDataTogether(i,:),'WriteMode','append');
    end

    filenamel = sprintf('%s_lo.h5','data');
    name = sprintf('/RecordB%g',i);
    if(i==1)
        hdf5write(filenamel,name,lowbandDataTogether(i,:));
    else
        hdf5write(filenamel,name,lowbandDataTogether(i,:),'WriteMode','append');
    end
    
end