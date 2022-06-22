function [phys] = TNC_ExtractNEV2IndChan(nevDS)

disp(' ');
disp('......................................................................');
disp(['Analyzing: ' nevDS.MetaTags.Filename]);

sUnInds  = find(nevDS.Data.Spikes.Unit>0);
numUnInds= numel(sUnInds);

eArray = double(nevDS.Data.Spikes.Electrode(sUnInds));
uArray = double(nevDS.Data.Spikes.Unit(sUnInds));
wArray = double(nevDS.Data.Spikes.Waveform(sUnInds,:));
tArray = double(nevDS.Data.Spikes.Timestamps(sUnInds));

phys.maxTime = max(tArray);

uCnt = 1;

for j = 1:128 % only check ephys channels at the moment. Only have a 128 channel system at the moment
    
    currEv = find(eArray==j);

    if numel(currEv)>0
        
        uArrayC = uArray(currEv);
        eArrayC = eArray(currEv);
        wArrayC = wArray(currEv,:);
        tArrayC = tArray(currEv);      
        
        uList   = unique(uArrayC);
        
        for k=1:numel(uList)
            currUsI = find(uArrayC==k);
            if numel(currUsI)>0
                disp(['     Assigning unit ' num2str(k) ' on electrode ' num2str(j) ' to the index ' num2str(uCnt)]);
                phys.unit(uCnt).el  = j;
                phys.unit(uCnt).id  = k;
                phys.unit(uCnt).ts  = tArrayC(currUsI);
                phys.unit(uCnt).wf  = wArrayC(currUsI,:);
                uCnt                = uCnt+1;
            end
        end
    end
    
end

phys.numUnits = uCnt-1;


disp('......................................................................');
disp(' ');
