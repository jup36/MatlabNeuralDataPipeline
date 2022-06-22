function Update_text_boxes(handles)
% Read and save the information introduced by the user in text boxes

USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Neurons = USER_DATA{2};
Electrode = USER_DATA{3};
Par_sim = USER_DATA{4};

% Manual single units firing rates
Neurons.rates_manual(1) = str2double(get(handles.rate1,'String'));
Neurons.rates_manual(2) = str2double(get(handles.rate2,'String'));
Neurons.rates_manual(3) = str2double(get(handles.rate3,'String'));
Neurons.rates_manual(4) = str2double(get(handles.rate4,'String'));
Neurons.rates_manual(5) = str2double(get(handles.rate5,'String'));

% Manual single units distance
Neurons.distance_manual(1) = str2double(get(handles.distance1,'String'));
if Neurons.distance_manual(1) < 0
    Neurons.distance_manual(1) = 0;
    set(handles.distance1,'String','0');
end
if Neurons.distance_manual(1) > 1
    Neurons.distance_manual(1) = 1;
    set(handles.distance1,'String','1');
end

Neurons.distance_manual(2) = str2double(get(handles.distance2,'String'));
if Neurons.distance_manual(2) < 0
    Neurons.distance_manual(2) = 0;
    set(handles.distance2,'String','0');
end
if Neurons.distance_manual(2) > 1
    Neurons.distance_manual(2) = 1;
    set(handles.distance2,'String','1');
end

Neurons.distance_manual(3) = str2double(get(handles.distance3,'String'));
if Neurons.distance_manual(3) < 0
    Neurons.distance_manual(3) = 0;
    set(handles.distance3,'String','0');
end
if Neurons.distance_manual(3) > 1
    Neurons.distance_manual(3) = 1;
    set(handles.distance3,'String','1');
end

Neurons.distance_manual(4) = str2double(get(handles.distance4,'String'));
if Neurons.distance_manual(4) < 0
    Neurons.distance_manual(4) = 0;
    set(handles.distance4,'String','0');
end
if Neurons.distance_manual(4) > 1
    Neurons.distance_manual(4) = 1;
    set(handles.distance4,'String','1');
end

Neurons.distance_manual(5) = str2double(get(handles.distance5,'String'));
if Neurons.distance_manual(5) < 0
    Neurons.distance_manual(5) = 0;
    set(handles.distance5,'String','0');
end
if Neurons.distance_manual(5) > 1
    Neurons.distance_manual(5) = 1;
    set(handles.distance5,'String','1');
end

% Manual single units selected
if Par_sim.manual == 1
    n(1) = get(handles.neuron1,'Value');
    n(2) = get(handles.neuron2,'Value');
    n(3) = get(handles.neuron3,'Value');
    n(4) = get(handles.neuron4,'Value');
    n(5) = get(handles.neuron5,'Value');
    Neurons.nmanual = n.*[1 2 3 4 5];
    Neurons.nmanual = Neurons.nmanual(Neurons.nmanual~=0);
end

USER_DATA{1} = Par_cube;  
USER_DATA{2} = Neurons;  
USER_DATA{3} = Electrode;  
USER_DATA{4} = Par_sim;

set(handles.neurocube_figure,'userdata',USER_DATA);