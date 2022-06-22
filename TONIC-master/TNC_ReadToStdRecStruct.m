function [tonicDataStructure] = TNC_ReadToStdRecStruct(dataStructure,unitArrayToLoad,CSid,USid,ELid,digFlag,extFlag)
% FUNCTION DETAILS: High level function that builds a single structure containing a standard set of analyses for every unit
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
%

lost = 1;
tonicDataStructure.events.US.ts     = 0;
tonicDataStructure.events.EL.ts     = 0;

if digFlag(1,1) == 0
    % analog events
    tonicDataStructure.events.CS.ts     = dataStructure.neurons{CSid}.timestamps.*1000;
else
    % digital events
    tonicDataStructure.events.CS.ts     = TNC_GetDigitalStampsNEX(dataStructure,CSid);
end

disp([ 'Number of trials: ' num2str(length(tonicDataStructure.events.CS.ts)) ]);

if extFlag==1

    tonicDataStructure.USexists = 0;
    tonicDataStructure.events.US.ts = 0;
    learningFlag                = 0;
    tmpTs                       = 0;

else
    
    tonicDataStructure.USexists = 1;
    learningFlag=1;

    if digFlag(1,2) == 0
        if USid<0
            tonicDataStructure.events.US.ts     = tonicDataStructure.events.CS.ts + (-1000.*USid);
        else
            tonicDataStructure.events.US.ts     = dataStructure.neurons{USid}.timestamps.*1000;
        end
    else
        tonicDataStructure.events.US.ts     = TNC_GetDigitalStampsNEX(dataStructure,USid);
    end
    
end

% Finding the licks and what not
if digFlag(1,3) == 0
    if ELid==0
        learning = 0;
        tmpTs    = 0;
    else
        tonicDataStructure.events.EL.ts     = dataStructure.neurons{ELid}.timestamps.*1000;
        tmpTs       = tonicDataStructure.events.EL.ts;
    end
else
    tonicDataStructure.events.EL.ts     = TNC_GetDigitalStampsNEX(dataStructure,ELid);
    if length(tonicDataStructure.events.EL.ts)>2
        tmpTs    = tonicDataStructure.events.EL.ts;
    else
        tmpTs    = 0;
    end
end
    
% created aligned response functions for the continuous every lick recording
numStamps   = length(tonicDataStructure.events.EL.ts);

if numStamps>2
    tonicDataStructure.events.EL.delta                  = zeros(round(tmpTs(numStamps,1)),1)';
    tonicDataStructure.events.EL.delta(round(tmpTs))    = 1;
else
    tonicDataStructure.events.EL.delta                  = zeros(1000,1)';    
end

disp(['CS: ' num2str(length(tonicDataStructure.events.CS.ts)) ' | US: ' num2str(length(tonicDataStructure.events.US.ts)) ' | EL: ' num2str(length(tonicDataStructure.events.EL.ts))]);

neuronCount = 1;

for i=1:size(unitArrayToLoad,2)
    
    ch = unitArrayToLoad(i);

    % store as ms spaced timestamp data
    tonicDataStructure.unit(neuronCount).ts     = dataStructure.neurons{ch}.timestamps.*1000;
    tonicDataStructure.unit(neuronCount).name   = dataStructure.neurons{ch}.name;
    
    % need to search through this and make sure that i have the correct corresponding waveform since there can be various forms (template or wf).
    for z=1:size(dataStructure.waves,1)
        if strcmp(dataStructure.waves{z}.name,[dataStructure.neurons{ch}.name '_wf']);
            tonicDataStructure.unit(neuronCount).wf     = dataStructure.waves{z}.waveforms;
            lost = 0; % found
        end
    end
    
    if lost
        for z=1:size(dataStructure.waves,1)
            if strcmp(dataStructure.waves{z}.name,[dataStructure.neurons{ch}.name '_template']);
                tonicDataStructure.unit(neuronCount).wf     = dataStructure.waves{z}.waveforms;
                lost = 0; % found
            end
        end
    end
    
    if lost
        tonicDataStructure.unit(neuronCount).wf = zeros(1,48);
        lost = 0; % found        
    end
    
%     % store some temporary data that makes the function easier to read
%     validTS     = tonicDataStructure.unit(neuronCount).ts>1;
%     tmpTs       = tonicDataStructure.unit(neuronCount).ts(validTS);
%     numStamps   = length(tmpTs);

    % calculate the properties of the isi distributions
    [isi] = TNC_QuantISI(tonicDataStructure.unit(neuronCount).ts);
    tonicDataStructure.unit(neuronCount).isi = isi;

    disp(['Processed unit ' num2str(unitArrayToLoad(i)) ' (' dataStructure.neurons{ch}.name ').....events: ' num2str(length(tonicDataStructure.unit(neuronCount).ts))]);
     
    neuronCount = neuronCount+1;
end
