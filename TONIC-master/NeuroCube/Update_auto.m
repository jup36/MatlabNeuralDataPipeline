function Update_auto(handles)

Update_text_boxes(handles)
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Neurons = USER_DATA{2};
Electrode = USER_DATA{3};
Par_sim = USER_DATA{4};

if Par_sim.manual == 1
    x_neurons = Neurons.coordinates(:,1)';
    y_neurons = Neurons.coordinates(:,2)';
    z_neurons = Neurons.coordinates(:,3)';
    distance = zeros(Electrode.nchannels,Neurons.nneurons);
    for i = 1 : Electrode.nchannels
        distance(i,:) = sqrt((x_neurons(1,:)-Electrode.coordinates(i,1)).^2+(y_neurons(1,:)-Electrode.coordinates(i,2)).^2+...
            (z_neurons(1,:)-Electrode.coordinates(i,3)).^2);
    end
    % Delete auto single units 
    if Electrode.nchannels == 1
        single_units_distance = (Electrode.diameter/2) + Par_sim.dsingle;
        single_units = ceil(find(distance<single_units_distance)./Electrode.nchannels);
        min_distance = (Electrode.diameter/2) + Par_cube.margin_electrode;
    else
        single_units_distance = abs(Electrode.coordinates(1,1)) + (Electrode.diameter/2) + Par_sim.dsingle;
        single_units = ceil(find(distance<single_units_distance)./Electrode.nchannels);
        min_distance = abs(Electrode.coordinates(1,1)) + (Electrode.diameter/2) + Par_cube.margin_electrode;
    end
    if ~isempty(single_units)
        distance(:,single_units) = [];
        x_neurons(single_units) = [];
        y_neurons(single_units) = [];
        z_neurons(single_units) = [];
        Neurons.id(single_units) = [];
        Neurons.Rates(single_units) = [];
        Neurons.spiketimes(single_units) = [];
    end
    % Include manual single units
    for i = 1 : length(Neurons.nmanual)
        single_aux = Neurons.nmanual(i);
        real_d = (single_units_distance-min_distance)*Neurons.distance_manual(single_aux) + min_distance;
        [x,y,z] = Manual_neuron(real_d);
        x_neurons = [x_neurons x];
        y_neurons = [y_neurons y];
        z_neurons = [z_neurons z];
        real_rate(i) = Neurons.rates_manual(single_aux);
    end
    Neurons.nneurons = length(x_neurons);
    Neurons.coordinates = [];
    Neurons.coordinates = [x_neurons' y_neurons' z_neurons'];
    
    if ~isempty(Neurons.nmanual)
        Models = randi(5,1,length(Neurons.nmanual));
        Parameters = randi(4,1,length(Neurons.nmanual));
        Neurons.id(end+1:end+length(Neurons.nmanual)) = 600 + (Models-1)*4 + Parameters;
        Neurons.Rates(end+1:end+length(Neurons.nmanual)) = real_rate;
    end
    
    tres = 1/(Par_sim.downsampling*Par_sim.sr);
    Refract = Par_sim.refract_period * 0.001 / tres;
    for i = 1 : length(Neurons.nmanual)
        index = length(Neurons.Rates)-length(Neurons.nmanual)+i;
        ISI = exprnd(Par_sim.sr*Par_sim.downsampling/Neurons.Rates(index),1,...
            round(Neurons.Rates(index)*Par_sim.duration));
        spiketimes = round(cumsum(ISI));
        CloseInds = diff(spiketimes) <= Refract; % Exclude potential overlapping spikes for a given neuron
        spiketimes(CloseInds) = [];
        Neurons.spiketimes{index} = spiketimes;
    end
end
USER_DATA{1} = Par_cube;  
USER_DATA{2} = Neurons;  
USER_DATA{3} = Electrode;  
USER_DATA{4} = Par_sim;

set(handles.neurocube_figure,'userdata',USER_DATA);
