function Update_plot(handles)
% Update the plot of the cube

USER_DATA = get(handles.neurocube_figure,'userdata');

Par_cube = USER_DATA{1};
Neurons = USER_DATA{2};
Electrode = USER_DATA{3};
x_neurons = Neurons.coordinates(:,1);
y_neurons = Neurons.coordinates(:,2);
z_neurons = Neurons.coordinates(:,3);
volume_close = find(Neurons.id >= 600);
volume_far = find(Neurons.id < 600);
distance_max = 1000*Par_cube.edge/2;
nneurons_far = length(volume_far);
nneurons_close = length(volume_close);
nneurons_far_plot = round(nneurons_far*Par_cube.prop_plot/100);
nneurons_close_plot = round(nneurons_close*Par_cube.prop_plot/100);

cla(handles.Cube_plot)

x_neurons_far = x_neurons(volume_far(1:nneurons_far_plot));
y_neurons_far = y_neurons(volume_far(1:nneurons_far_plot));
z_neurons_far = z_neurons(volume_far(1:nneurons_far_plot));
plot3(handles.Cube_plot,x_neurons_far,y_neurons_far,z_neurons_far,'^c','markersize',3,'linewidth',3)
set(handles.Cube_plot,'XGrid','on','YGrid','on','ZGrid','on')
hold(handles.Cube_plot,'on')

x_neurons_close = x_neurons(volume_close(1:nneurons_close_plot));
y_neurons_close = y_neurons(volume_close(1:nneurons_close_plot));
z_neurons_close = z_neurons(volume_close(1:nneurons_close_plot));
plot3(handles.Cube_plot,x_neurons_close,y_neurons_close,z_neurons_close,'^r','markersize',3,'linewidth',3)
set(handles.Cube_plot,'XGrid','on','YGrid','on','ZGrid','on')
axis(handles.Cube_plot,[-distance_max distance_max -distance_max distance_max -distance_max distance_max])
hold(handles.Cube_plot,'on')

set(handles.neurocube_figure,'CurrentAxes',handles.Cube_plot)
for i = 1 : Electrode.nchannels
    h = filledCircle(Electrode.coordinates(i,:),Electrode.diameter/2,100,'k');
end