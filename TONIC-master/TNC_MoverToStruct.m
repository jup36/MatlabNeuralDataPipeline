function [tonicDataStructure] = TNC_MoverToStruct(dataStructure,unitArrayToLoad)
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

tonicDataStructure.events.TRIAL.ts     = TNC_GetDigitalStampsNEX(dataStructure,'-0006');
tonicDataStructure.events.REWARD.ts       = TNC_GetDigitalStampsNEX(dataStructure,'-0007');
tonicDataStructure.events.LICKS.ts      = TNC_GetDigitalStampsNEX(dataStructure,'-0008');

disp([ 'Number of trials: '   num2str(length(tonicDataStructure.events.TRIAL.ts))  ]);
disp([ 'Number of licks: '    num2str(length(tonicDataStructure.events.LICKS.ts)) ]);

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

    % calculate the properties of the isi distributions
    [isi] = TNC_QuantISI(tonicDataStructure.unit(neuronCount).ts);
    tonicDataStructure.unit(neuronCount).isi = isi;

    disp(['Processed unit ' num2str(unitArrayToLoad(i)) ' (' dataStructure.neurons{ch}.name ').....events: ' num2str(length(tonicDataStructure.unit(neuronCount).ts))]);
     
    neuronCount = neuronCount+1;
end
