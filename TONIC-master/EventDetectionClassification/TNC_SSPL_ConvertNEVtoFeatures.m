function [] = TNC_SSPL_ConvertNEVtoFeatures(fileNameStr,chunk,arrayType,targetName)

%% PARAMETERS
    upRatio = 5;
    debug = 0;
    
%% LOAD the file data

    [nevData]   = TNC_LoadData(0, 0, fileNameStr);
    numSegs     = ceil(max(nevData.Data.Spikes.Timestamps) ./ (chunk.*30000) );
    allElec     = unique(nevData.Data.Spikes.Electrode);
    trodeList   = allElec(find(allElec>0 & allElec<=128));
    numTrodes   = numel(trodeList);
    numSamps    = numel(nevData.Data.Spikes.Waveform(1,:));

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

    disp('_________________________________________________________________ ');
    disp(' ');
    disp(['Analyzing events over the range of ' num2str((k-1)*chunk) ' to ' num2str(k*chunk) ' seconds.']);
    disp('_________________________________________________________________ ');

    for i=1:max(trodeList)
         
        clear spikes

        %________________________________________________________________________
        %________________________________________________________________________
        % GRAB SPIKES
                    
            spikes.inds = find(nevData.Data.Spikes.Electrode == i & nevData.Data.Spikes.Timestamps>begT & nevData.Data.Spikes.Timestamps<=endT);
            
            if numel(spikes.inds)>4

                spikes.wfs(numel(spikes.inds)).values = zeros(1,1,'int16');
                for zz=1:numel(spikes.inds)
                    spikes.wfs(zz).values  = int16(nevData.Data.Spikes.Waveform(spikes.inds(zz),:)); 
                end
                spikes.ts   = int16(nevData.Data.Spikes.Timestamps(spikes.inds,:));       
                spikes.num  = numel(spikes.inds);
                spikes.winL = 15;
                spikes.winR = 32;

                sessionStruct.seg(k).shank(i).wfs   = spikes.wfs;
                sessionStruct.seg(k).shank(i).inds  = spikes.inds;

        %________________________________________________________________________
        %________________________________________________________________________
        % CALC SCALAR QUANT FOR SPIKES
        
                disp(['Processing spikes on elec' num2str(i) '...']);          

                [features]   = TNC_SSPL_EventQuantSE(spikes,[]);            

                if debug
                    figure(200); plotmatrix(features.params(:,12:15)); pause(0.1);
                end
                
                featStruct.seg(k).shank(i).ts     = round(spikes.ts./30);
                featStruct.seg(k).shank(i).id     = zeros(numel(spikes.inds),1);
                featStruct.seg(k).shank(i).inds   = spikes.inds;
                featStruct.seg(k).shank(i).params = features.params;
                featStruct.paramNames             = features.paramNames;

            else
                
                sessionStruct.seg(k).shank(i).wfs(1).values = zeros(1,numSamps);
                sessionStruct.seg(k).shank(i).inds  = 0;
                featStruct.seg(k).shank(i).id     = 0;
                featStruct.seg(k).shank(i).inds   = 0;
                featStruct.seg(k).shank(i).ts     = 0;
                featStruct.seg(k).shank(i).params = 0;
                
            end
    end
    
end

%% SAVE FEATURE STRUCTURE

    offset                  = 0;
    featStruct.chunk        = chunk;
    featStruct.arrayType    = arrayType;
    disp(' ');
    disp('_________________________________________________________________ ');
    disp(['save ' targetName '_ft featStruct']);
    save([targetName '_ft.mat'],'featStruct');
    
%% SAVE SESSION STRUCTURE
 
disp(['save ' targetName '_ss sessionStruct']);
save([targetName '_ss.mat'],'sessionStruct');