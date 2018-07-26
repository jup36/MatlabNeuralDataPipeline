function [DDS] = TNC_ConvertSortedNEVtoDDS(fileNameStr,outResolution,sessNum)
% FUNCTION DETAILS: Convert an NEV file with sorted spikes to the DudmanLabDataStructure (DDS)
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: josh@dudmanlab.org
% CONTRIBUTIONS: www.dudmanlab.org
% _________________________________________________________________________
%
% INPUT PARAMETERS:
%   fileNameStr - path to .nev file to be loaded
%   outResolution - in kHz (default is 1 kHz - ms resolution)
%   sessNum - datastructures are organized into sessions at the top level
% 
% 

j=1;

scaleFactor = 30 ./ outResolution;

EV_Active = TNC_LoadData(0, 0, fileNameStr);

avs_ind     = find(EV_Active.Data.Spikes.Unit>0 & EV_Active.Data.Spikes.Unit<100);
avs_trode   = EV_Active.Data.Spikes.Electrode(avs_ind);
avs_allTrode= unique(avs_trode);

for k=1:numel(avs_allTrode)
    
    currTrode = avs_allTrode(k);
    
    if currTrode<128
        
        trode_inds = find(EV_Active.Data.Spikes.Electrode==currTrode);
        unit_nums = unique(EV_Active.Data.Spikes.Unit(trode_inds));
        
        for m=1:numel(unit_nums)

            curr_unit = unit_nums(m);
            
            if curr_unit>0 & curr_unit<20

                acs_inds                            = find(EV_Active.Data.Spikes.Unit==curr_unit & EV_Active.Data.Spikes.Electrode==currTrode);
                DDS.session(sessNum).unit(j).rts    = EV_Active.Data.Spikes.Timestamps(acs_inds);
                DDS.session(sessNum).unit(j).ts     = round( DDS.session(sessNum).unit(j).rts ./ scaleFactor );
                DDS.session(sessNum).unit(j).el     = currTrode;                
                DDS.session(sessNum).unit(j).wf     = EV_Active.Data.Spikes.Waveform(acs_inds,:)';
                DDS.session(sessNum).unit(j).ui     = curr_unit;
                DDS.session(sessNum).unit(j).da     = 0; % to be set later...

                disp([num2str(j) ' >>> elec:' num2str(currTrode) ' | unit:' num2str(curr_unit) ' ... ' num2str(numel(DDS.session(sessNum).unit(j).ts)) ' spks' ]);

                j = j+1;
                
            end                    
        end                
    end
end

DDS.session(sessNum).num_units = j-1;

disp(' ');
disp('_____________________________________________________________________________________');
disp('_____________________________________________________________________________________');
disp(' ');
disp(['Extracted ' num2str(j-1) ' sorted units from the file: ' fileNameStr]);
disp('_____________________________________________________________________________________');
disp('_____________________________________________________________________________________');
disp(' ');
