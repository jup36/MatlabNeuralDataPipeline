%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function TNC_MM_UpdateNeuroCube(varargin)      

    debug = 1;

    if nargin ==  3
        ncFile = varargin{1};
        Nunits = varargin{2};
        dlim   = varargin{3};
    else
        disp('Usage: TNC_MM_UpdateNeuroCube(Neurocube_file, NumUnits, dlim)');
        return;
    end
    ncData = load(ncFile); 
    Neurons.nneurons = ncData.Neurons.nneurons;
    Neurons.coordinates = ncData.Neurons.coordinates;
    Neurons.nmanual = [];
    Neurons.distance_manual = ncData.Neurons.distance_manual;
    Neurons.rates_manual    = ncData.Neurons.rates_manual;
    Neurons.id              = ncData.Neurons.id;
    Neurons.Rates           = ncData.Neurons.Rates;
    Neurons.spiketimes      = ncData.Neurons.spiketimes;
    Neurons_aux = ncData.Neurons_aux;
    n = 1;
    for i=1:length(ncData.Neurons.spiketimes)
%       disp(['i=' num2str(i) ' length(ncData.Neurons.spiketimes(i)=' length(ncData.Neurons.spiketimes(i))]);
%       if n <= length(XC) && length(ncData.Neurons.spiketimes{i}) > 0
        dist = sqrt(Neurons.coordinates(i,1)^2 + Neurons.coordinates(i,2)^2 + Neurons.coordinates(i,3)^2);
        if length(Neurons.spiketimes{i}) > 0 && dist < dlim && n <= Nunits
            n = n + 1;
            if debug
                disp(['i=' num2str(i) ' length(spiketimes(i))=' num2str(length(Neurons.spiketimes{i}))...
                      ' dist=' num2str(dist) ' X,Y,Z=' num2str(Neurons.coordinates(i,1)) ' ' num2str(Neurons.coordinates(i,2)) ' ' num2str(Neurons.coordinates(i,3)) ]);
            end
        else
            Neurons.spiketimes{i} = []; 
        end 
    end
    outputName  = strcat('NeuroCubeUpdated_', num2str(Nunits));
    disp(['save ' outputName  ' Neurons Neurons_aux ']);
    eval(['save ' outputName  ' Neurons Neurons_aux ']);


