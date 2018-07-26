function [] = TNC_SSPL_ConvertNEVtoFeatures(fileNameStr,chunk,arrayType,targetName)

%% PARAMETERS
    upRatio = 5;
    debug = 0;
    
%% LOAD the file data

    [nevData] = TNC_LoadData(0, 0, fileNameStr);
    numSegs = floor(max(nevData.Data.Spikes.Timestamps) ./ (chunk.*30000) );
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

    disp('_________________________________________________________________ ');
    disp(' ');
    disp(['Analyzing events over the range of ' num2str((k-1)*chunk) ' to ' num2str(k*chunk) ' seconds.']);
    disp('_________________________________________________________________ ');

    for i=1:numTrodes
                
        %________________________________________________________________________
        %________________________________________________________________________
        % GRAB SPIKES
                    
            spikes.inds = find(nevData.Data.Spikes.Electrode == trodeList(i) & nevData.Data.Spikes.Timestamps>begT & nevData.Data.Spikes.Timestamps<=endT);
            
            if numel(spikes.inds)>4
                spikes.wfs  = double(nevData.Data.Spikes.Waveform(spikes.inds,:));        
                spikes.ts   = double(nevData.Data.Spikes.Timestamps(spikes.inds,:));       
                spikes.num  = numel(spikes.inds);
                spikes.winL = 15;
                spikes.winR = 32;

                sessionStruct.seg(k).shank(i).wfs   = spikes.wfs;
                sessionStruct.seg(k).shank(i).inds  = spikes.inds;

%             disp(['Align spikes from elec ' num2str(trodeList(i)) '.']);
%             [events]    = TNC_EventAlign(spikes,upRatio,0);

        %________________________________________________________________________
        %________________________________________________________________________
        % CALC SCALAR QUANT FOR SPIKES
        
%                 disp(['Quantify spikes from elec ' num2str(trodeList(i)) '.']);          

                [features]   = TNC_EventQuantSE(spikes,[]);            

                if debug
                    figure(200); plotmatrix(features.params(:,12:15)); pause(0.1);
                end
                
                featStruct.seg(k).shank(i).inds   = spikes.inds;
                featStruct.seg(k).shank(i).ts     = round(spikes.inds./30);
                featStruct.seg(k).shank(i).id     = zeros(numel(spikes.inds),1);
                featStruct.seg(k).shank(i).params = features.params;
                featStruct.paramNames             = features.paramNames;

            else
                
                sessionStruct.seg(k).shank(i).wfs   = [];
                sessionStruct.seg(k).shank(i).inds  = [];
                featStruct.seg(k).shank(i).inds   = [];
                featStruct.seg(k).shank(i).ts     = [];
                featStruct.seg(k).shank(i).params = [];
                
            end
    end
    
end

%% SAVE FEATURE STRUCTURE

    offset                  = 0;
    featStruct.chunk        = chunk;
    featStruct.arrayType    = arrayType;

    disp(['save ' targetName '_ft featStruct']);
    save([targetName '_ft.mat'],'featStruct');
    
%% SAVE SESSION STRUCTURE
 
disp(['save ' targetName '_ss sessionStruct']);
save([targetName '_ss.mat'],'sessionStruct');