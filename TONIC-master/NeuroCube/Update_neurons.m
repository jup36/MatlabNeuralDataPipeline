function Update_neurons(handles)
% Update the classification between close and distant neurons

USER_DATA = get(handles.neurocube_figure,'userdata');

Par_cube = USER_DATA{1};
Neurons = USER_DATA{2};
Electrode = USER_DATA{3};
Par_sim = USER_DATA{4};
x_neurons = Neurons.coordinates(:,1)';
y_neurons = Neurons.coordinates(:,2)';
z_neurons = Neurons.coordinates(:,3)';
volume_close_ini = find(Neurons.id >= 600);
volume_far = find(Neurons.id < 600);

distance = zeros(Electrode.nchannels,Neurons.nneurons);
for i = 1 : Electrode.nchannels
    distance(i,:) = sqrt((x_neurons(1,:)-Electrode.coordinates(i,1)).^2+(y_neurons(1,:)-Electrode.coordinates(i,2)).^2+...
        (z_neurons(1,:)-Electrode.coordinates(i,3)).^2);
end

volume_close = ceil(find(distance<=(Par_sim.dlim+Electrode.diameter/2))./Electrode.nchannels);   % Neurons to be simulated
repeated = diff(volume_close)==0;
volume_close(repeated) = [];

neurons_disappear=volume_close_ini(~ismember(volume_close_ini,volume_close));
neurons_appear=volume_close(~ismember(volume_close,volume_close_ini));

% Neurons that were close and now are distant

if ~isempty(neurons_disappear)
    load spike_shapes.mat
    SpikeShapes = NormShapes./repmat(max(NormShapes,[],2),1,size(NormShapes,2));
    clear Shapes
    [nshapes,~] = size(SpikeShapes);
    Neurons.id(neurons_disappear) = randi(nshapes,1,length(neurons_disappear));
%     Neurons.spiketimes{neurons_disappear} = [];
end

% Neurons that were distant and now are close

if ~isempty(neurons_appear)
    Models = randi(5,1,length(neurons_appear));
    Parameters = randi(4,1,length(neurons_appear));
    Neurons.id(neurons_appear) = 600 + (Models-1)*4 + Parameters;
    
    downsampling = 4;
    tres = 1/(downsampling*Par_sim.sr);
    Ns1 = round(Par_sim.duration / tres);
    Refract = Par_sim.refract_period * 0.001 / tres;
    for i = 1 : length(neurons_appear)
        ISI = exprnd(Par_sim.sr*downsampling/Neurons.Rates(neurons_appear(i)),1,round(Neurons.Rates(neurons_appear(i))*Par_sim.duration));
        spiketimes = round(cumsum(ISI));
%         spiketimes = find(binornd(1,Neurons.Rates(neurons_appear(i))*tres,1,Ns1));
        CloseInds = diff(spiketimes) <= Refract; % Exclude potential overlapping spikes for a given neuron
        spiketimes(CloseInds) = [];
        Neurons.spiketimes{neurons_appear(i)} = spiketimes;
    end
end

% Neurons that are too close to the electrode and must be deleted

min_distance = (Electrode.diameter/2) + Par_cube.margin_electrode; 
too_close = ceil(find(distance<min_distance)./Electrode.nchannels);
repeated = diff(too_close)==0;
too_close(repeated) = [];
if ~isempty(too_close)
    Neurons.coordinates(too_close,:) = [];
    Neurons.id(too_close) = [];
    Neurons.Rates(too_close) = [];
    Neurons.spiketimes(too_close) = [];
    Neurons.nneurons = Neurons.nneurons - length(too_close);
end

USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);