function [] = TNC_SSP_ConvertNEVtoFeatures(fileNameStr,chunk)

%% PARAMETERS
    upRatio = 5;

%% LOAD the file data

    [nevData] = TNC_LoadData(0, 0, fileNameStr);
    numSegs = floor(max(nevData.Data.Spikes.Timestamps) ./ chunk);
    allElec = unique(nevData.Data.Spikes.Electrode);
    trodeList = allElec(find(allElec>0 & allElec<=128));
    numTrodes = numel(trodeList);

    disp(' ');
    disp(' ');
    disp(['Data will be loaded and processed as ' num2str(numSegs) ' x ' num2str(chunk) 'sec segments.']);
    disp(['...Found ' num2str(numTrodes) ' electrodes.']);
    disp(' ');
    disp(' ');
    
%% COMPUTE features from waveforms chunk by chunk

for k=1:numSegs
    
    begT = (k-1)*chunk*30000;
    endT = k*chunk*30000;
    
    disp(['Analyzing events over the range of ' num2str((k-1)*chunk) ' to ' num2str(k*chunk) ' seconds.']);

    for i=1:numTrodes
                
        %________________________________________________________________________
        %________________________________________________________________________
        % ALIGN SPIKES
        
            disp(['Align spikes from elec ' num2str(trodeList(i)) '.']);
            spikes.inds = find(nevData.Data.Spikes.Electrode == trodeList(i) & nevData.Data.Spikes.Timestamps>begT & nevData.Data.Spikes.Timestamps<=endT);
            spikes.wfs  = nevData.Data.Spikes.Waveform(events.inds,:);        
            spikes.ts   = nevData.Data.Spikes.Timestamps(events.inds,:);       
            spikes.num  = numel(events.inds);
            spikes.winL = 15;
            spikes.winR = 32;

            [events]    = TNC_EventAlign(spikes,upRatio,0);

        %________________________________________________________________________
        %________________________________________________________________________
        % CALC SCALAR QUANT FOR SPIKES
        
            disp(['Quantify spikes from elec ' num2str(trodeList(i)) '.']);
        
            [quant]   = TNC_EventQuant(events,'scalar','dot',5);
            tmpParams = [quant.scalar.scl.values(:,[1,2,3,8,9])];

            [quant]   = TNC_EventQuant(events,'pca','dot',5);
            tmpParams = [tmpParams , quant.pca.dot.values];            

            featStruct.seg(k).shank(i).inds   = events.inds;
            featStruct.seg(k).shank(i).ts     = round(events.inds./30);
            featStruct.seg(k).shank(i).params = tmpParams;
            
    end
    
end

%% Save the output file

    featStruct.paramnames = {'min','max','nEng','eng','int','pc1','pc2','pc3','pc4','pc5'};

    disp(['save ' fileName(1:totChar-4) '_ft featStruct']); 
    eval(['save ' fileName(1:totChar-4) '_ft featStruct']); 
    