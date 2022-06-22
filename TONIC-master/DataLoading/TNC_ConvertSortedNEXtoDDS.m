function [DDS] = TNC_ConvertSortedNEXtoDDS(fileNameStr,sessNum)
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
%   fileNameStr - path to .nex file to be loaded
%   outResolution - in kHz (default is 1 kHz - ms resolution)
%   sessNum - datastructures are organized into sessions at the top level
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD SORTED SPIKE DATA FROM NEUROEXPLORER     
    RD_Active = TNC_LoadData(0, 0, fileNameStr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

    propUA = []; propCS = 0; propEL = 0; clear propNames; p = 1;

    for m = 1:size(RD_Active.waves,1)    
        
        if numel(findstr(RD_Active.waves{m}.name,'Chan'))>0 % electrode
            newtype=0;

            testi = findstr(RD_Active.waves{m}.name,'u');
            if size(testi,1) == 0
                propUA = [propUA,m];
                propNames.names(p).str = RD_Active.waves{m}.name;
                disp(['Found unit ' num2str(p) ': ' RD_Active.waves{m}.name]);
                p = p+1;
            end
            
        elseif numel(findstr(RD_Active.waves{m}.name,'elec'))>0

            newtype=1;

%             if numel(findstr(RD_Active.waves{m}.name,'_ts'))>0 && numel(findstr(RD_Active.waves{m}.name,'_u'))>0
            tmpInd = findstr(RD_Active.waves{m}.name,'_u');
            if numel(tmpInd)>0 && str2num(RD_Active.waves{m}.name(tmpInd+2))~=0
                propUA = [propUA,m];
                propNames.names(p).str = RD_Active.waves{m}.name;
                disp(['Found unit ' num2str(p) ': ' RD_Active.waves{m}.name]);
                p = p+1;                
            end
        end
    end

    p=p-1;

    if newtype
        % Write into the pop data structure
        for j=1:p
            i = propUA(j);
            DDS.session(sessNum).unit(j).ts = round(RD_Active.waves{i}.timestamps.*1000);
            DDS.session(sessNum).unit(j).el = str2double(RD_Active.waves{i}.name(5:6));        
            DDS.session(sessNum).unit(j).wf = RD_Active.waves{i}.waveforms;
            DDS.session(sessNum).unit(j).ui = str2double( RD_Active.waves{i}.name(findstr(RD_Active.waves{i}.name,'_u')+2) );
            DDS.session(sessNum).unit(j).da = 0; % to be set later...

            [isi] = TNC_QuantISI(DDS.session(sessNum).unit(j).ts);
            DDS.session(sessNum).unit(j).isi = isi; 

            disp(['Unit ' num2str(j) ': elec' num2str(DDS.session(sessNum).unit(j).el) '_unit' num2str(DDS.session(sessNum).unit(j).ui)]);
        end
    else
        % Write into the pop data structure
        for j=1:p
            i = propUA(j);
            DDS.session(sessNum).unit(j).ts = round(RD_Active.waves{i}.timestamps.*1000);
            if str2num(propNames.names(j).str(6))==0
                DDS.session(sessNum).unit(j).el = str2num(propNames.names(j).str(7));        
            else
                DDS.session(sessNum).unit(j).el = str2num(propNames.names(j).str(6:7));        
            end
            DDS.session(sessNum).unit(j).wf = RD_Active.waves{i}.waveforms;
            DDS.session(sessNum).unit(j).ui = RD_Active.waves{i}.unitNumber;
            DDS.session(sessNum).unit(j).da = 0; % to be set later...

            [isi] = TNC_QuantISI(DDS.session(sessNum).unit(j).ts);
            DDS.session(sessNum).unit(j).isi = isi; 

            disp(propNames.names(j).str);
        end
        
    end
    
    DDS.session(sessNum).num_units = p;

disp(' ');
disp('_____________________________________________________________________________________');
disp('_____________________________________________________________________________________');
disp(' ');
disp(['Extracted ' num2str(p) ' sorted units from the file: ' fileNameStr]);
disp('_____________________________________________________________________________________');
disp('_____________________________________________________________________________________');
disp(' ');
    