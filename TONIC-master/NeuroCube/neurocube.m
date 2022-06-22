function varargout = neurocube(varargin)
% neurocube M-file for neurocube.fig
%      neurocube, by itself, creates a new neurocube or raises the existing
%      singleton*.
%
%      H = neurocube returns the handle to a new neurocube or the handle to
%      the existing singleton*.
%
%      neurocube('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in neurocube.M with the given input arguments.
%
%      neurocube('Property','Value',...) creates a new neurocube or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before neurocube_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to neurocube_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help neurocube

% Last Modified by GUIDE v2.5 02-Jan-2013 16:02:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @neurocube_OpeningFcn, ...
                   'gui_OutputFcn',  @neurocube_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before neurocube is made visible.
function neurocube_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to neurocube (see VARARGIN)

% Choose default command line output for neurocube
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes neurocube wait for user response (see UIRESUME)
% uiwait(handles.neurocube_figure);


% --- Outputs from this function are returned to the command line.
function varargout = neurocube_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

addpaths;                           % Add the sub-directories to the path

Par_cube.edge = 1;                  % Edge of the cube in mm
%Par_cube.density = 300000;          % Neurons/mm^3 in the cube
Par_cube.density = 50000;         
Par_cube.prop_firing = 7;           % Proportion of firing neurons in %
Par_cube.firing_rates = 1;          % 1 = exponential; 2 = uniform
Par_cube.firing_rates_params = ...
    [0.21 1 0 0.5 5];               % Exponential: k, sigma, theta; Uniform: min_rate, max_rate
Par_cube.margin_electrode = 10;     % The neurons inside a sphere with radius = radius_electrode + margin_electrode are deleted
Par_cube.prop_plot = 1;             % Proportion of neurons (in %) plotted

Neurons.nneurons = round(Par_cube.density*(Par_cube.edge^3)*(0.01*Par_cube.prop_firing));
                                    % Number of neurons to be generated in the cube
Neurons.coordinates = [];           % Coordinates of all the neurons in the cube
Neurons.nmanual = [];               % List of selected manual single units
Neurons.distance_manual = ...
    [0.5 0.5 0.5 0.5 0.5];          % Distances between the electrode and the manual single units
Neurons.rates_manual = [5 5 5 5 5]; % Firing rates of the manual single units

Electrode.nchannels = 1;            % Number of electrodes
Electrode.diameter = 20;            % Diameter of the electrodes in microns
Electrode.coordinates = [0 0 0];    % Coordinate of the single channel electrode
Electrode.separation = 50;          % Separation between electrodes in a tetrode

Electrode.tetrode = ...
    [1 1 0;1 -1 0;-1 1 0;-1 -1 0];  % Coordinates of the tetrode in microns
%   [0 0 0;1 0  0; 0 1 0; 1  1 0];  % Coordinates of the tetrode in microns

Par_sim.duration = 30;              % Duration of the simulation in seconds
Par_sim.sr = 24000;                 % Sampling rate
Par_sim.downsampling = 4;
Par_sim.refract_period = 2;         % Refractory period (in milliseconds)
Par_sim.dlim = 150;                 % Distance from the electrode to the limit of zone A (in microns)
Par_sim.dsingle = 100;              % Distance from the electrode that defines the volume of manually defined single units
Par_sim.manual = 0;                 % Flag that indicates that the single units are manually pre-defined (manual = 1)
Par_sim.std_th_noise = 1;           % Standard deviation of thermal noise in microvolts
Par_sim.grid_npoint = 25;           % Number of points in the grid used to calculate spike shapes across the electrode surface

Par_sim.point_source      = 1;
Par_sim.five_models       = 0;
Par_sim.orig_code         = 0;
Par_sim.keep_volume_far   = 0;
Par_sim.five_parameters   = 0;
Par_sim.orig_firing_rates = 0;
Par_sim.filter_data       = 0;

all_neuron_models         = {'d151' 'pyramidal264L' '51-2a' '51-2b' 're80'};
neuron_model_id           = 5;

if Par_sim.orig_code || Par_sim.five_models
    Par_sim.neuron_models = all_neuron_models';
else
    i = neuron_model_id;
    Par_sim.neuron_models = {all_neuron_models{i};all_neuron_models{i};...
                             all_neuron_models{i};all_neuron_models{i};...
                             all_neuron_models{i}};
end


USER_DATA = get(handles.neurocube_figure,'userdata');
USER_DATA{1} = Par_cube;            % Parameters of the cube structure
USER_DATA{2} = Neurons;             % Neurons structure
USER_DATA{3} = Electrode;           % Electrode structure
USER_DATA{4} = Par_sim;             % Simulation parameters structure

set(handles.neurocube_figure,'userdata',USER_DATA)

guidata(hObject, handles);


% --- Executes on button press in Generate.
function Generate_Callback(hObject, eventdata, handles)
% hObject    handle to Generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Update_text_boxes(handles)
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Neurons = USER_DATA{2};
Electrode = USER_DATA{3};
Par_sim = USER_DATA{4};

set(handles.Progress_text,'String','Generating cube...');
set(handles.Progress_text,'Visible','On');
pause(0.1)


% Calculate random neuron positions
% GD
disp(['Electrode.nchannels=' num2str(Electrode.nchannels)]);

keep_volume_close = 1;
distance_max = 1000*Par_cube.edge/2;    % Distance in microns from the center to the cube to the edge
if Par_sim.orig_code
    Neurons.nneurons = round(Par_cube.density*(Par_cube.edge^3)*(0.01*Par_cube.prop_firing));
    Neurons.x_neurons = randi([-distance_max distance_max],1,Neurons.nneurons);   % x coordinate of neurons in microns
    Neurons.y_neurons = randi([-distance_max distance_max],1,Neurons.nneurons);   % y coordinate of neurons in microns
    Neurons.z_neurons = randi([-distance_max distance_max],1,Neurons.nneurons);   % z coordinate of neurons in microns
else
    Neurons.nneurons_far = 0;
    if Par_sim.keep_volume_far
        % Generate random neurons
        Neurons.nneurons = round(Par_cube.density*(Par_cube.edge^3)*(0.01*Par_cube.prop_firing));
        disp(['Neurons.nneurons=' num2str(Neurons.nneurons)]);
        Neurons.x_neurons = randi([-distance_max distance_max],1,Neurons.nneurons);   % x coordinate of neurons in microns
        Neurons.y_neurons = randi([-distance_max distance_max],1,Neurons.nneurons);   % y coordinate of neurons in microns
        Neurons.z_neurons = randi([-distance_max distance_max],1,Neurons.nneurons);   % z coordinate of neurons in microns
        distance0 = zeros(Electrode.nchannels,Neurons.nneurons);
        for i = 1 : Electrode.nchannels
            distance0(i,:) = sqrt((Neurons.x_neurons(1,:)-Electrode.coordinates(i,1)).^2+...
                         (Neurons.y_neurons(1,:)-Electrode.coordinates(i,2)).^2+...
                         (Neurons.z_neurons(1,:)-Electrode.coordinates(i,3)).^2);
        end
        
        % Eliminate neurons in close volume
        volume_close = ceil(find(distance0<=(Par_sim.dlim+Electrode.diameter/2))./Electrode.nchannels);
        disp(['size(volume_close)=' num2str(size(volume_close))]);
        disp(['size(Neurons.x_neurons)=' num2str(size(Neurons.x_neurons))]);
        Neurons.x_neurons(volume_close) = [];
        disp(['size(Neurons.x_neurons_far)=' num2str(size(Neurons.x_neurons))]);
        Neurons.y_neurons(volume_close) = [];
        disp(['size(Neurons.y_neurons_far)=' num2str(size(Neurons.y_neurons))]);
        Neurons.z_neurons(volume_close) = [];
        disp(['size(Neurons.z_neurons_far)=' num2str(size(Neurons.z_neurons))]);
        Neurons.nneurons_far = size(Neurons.x_neurons,2); 
        disp(['Neurons.nneurons_far=' num2str(Neurons.nneurons_far)]);
    else
        Neurons.x_neurons = [];
        Neurons.y_neurons = [];
        Neurons.z_neurons = []
    end
    
    % Read input data from text file 'NeuroCubeInput.txt'
    if keep_volume_close
        i = 1;
        fid = fopen('NeuroCubeInput.txt','r'); %# open csv file for reading
        while ~feof(fid)
            line = fgets(fid); %# read line by line
            coords(i,1:3) = sscanf(line,'%f %f %f'); %# sscanf can read only numeric data :(
            i = i+1;
        end
        fclose(fid);
        disp(['Neurons.nneurons_close=' num2str(size(coords,1))]);
        Neurons.x_neurons = [ Neurons.x_neurons coords(:,1)' ];
        Neurons.y_neurons = [ Neurons.y_neurons coords(:,2)' ];
        Neurons.z_neurons = [ Neurons.z_neurons coords(:,3)' ];
        Neurons.nneurons = size(Neurons.x_neurons, 2);
    end
end
disp(' ');
disp(['Neurons.nneurons_total=' num2str(Neurons.nneurons)]);
disp(['size(Neurons.x_neurons)=' num2str(size(Neurons.x_neurons))]);
distance = zeros(Electrode.nchannels,Neurons.nneurons);
for i = 1 : Electrode.nchannels
    distance(i,:) = sqrt((Neurons.x_neurons(1,:)-Electrode.coordinates(i,1)).^2+...
                         (Neurons.y_neurons(1,:)-Electrode.coordinates(i,2)).^2+...
                         (Neurons.z_neurons(1,:)-Electrode.coordinates(i,3)).^2);
end

% GD
if Par_sim.orig_code
    if Neurons.nneurons == 0
      nlab = imread('filelist.xlj','jpg'); 
      scrsz = get(0,'ScreenSize'); figure('color','k','Position',[scrsz(1)+scrsz(3)/4 scrsz(2)+scrsz(4)/4.5 scrsz(3)/2 scrsz(4)/1.5]); 
      image(nlab); axis equal; axis image; axis off; set(gcf,'NumberTitle','off');
    end
end

% Delete neurons damaged by the electrode

min_distance = (Electrode.diameter/2) + Par_cube.margin_electrode; 
too_close = ceil(find(distance<min_distance)./Electrode.nchannels);
%disp('too close=');
%too_close
if Par_sim.orig_code
    distance(:,too_close) = [];
    Neurons.x_neurons(too_close) = [];
    Neurons.y_neurons(too_close) = [];
    Neurons.z_neurons(too_close) = [];
end

% If manual mode selected, delete auto single units and include manual ones

if Par_sim.manual == 1
    distance_auto = distance;
    x_neurons_auto = Neurons.x_neurons;
    y_neurons_auto = Neurons.y_neurons;
    z_neurons_auto = Neurons.z_neurons;
    if Electrode.nchannels == 1
        single_units_distance = (Electrode.diameter/2) + Par_sim.dsingle;
        single_units = ceil(find(distance<single_units_distance)./Electrode.nchannels);
        min_distance = (Electrode.diameter/2) + Par_cube.margin_electrode;
    else
        single_units_distance = abs(Electrode.coordinates(1,1)) + (Electrode.diameter/2) + Par_sim.dsingle;
        single_units = ceil(find(distance<single_units_distance)./Electrode.nchannels);
        min_distance = abs(Electrode.coordinates(1,1)) + (Electrode.diameter/2) + Par_cube.margin_electrode;
    end
    disp(['single_units=' num2str(single_units)]);
    if Par_sim.orig_code
        distance(:,single_units) = [];
        Neurons.x_neurons(single_units) = [];
        Neurons.y_neurons(single_units) = [];
        Neurons.z_neurons(single_units) = [];
    end
    for i = 1 : length(Neurons.nmanual)
        single_aux = Neurons.nmanual(i);
        real_d = (single_units_distance-min_distance)*Neurons.distance_manual(single_aux) + min_distance;
        [x,y,z] = Manual_neuron(real_d);
        Neurons.x_neurons = [Neurons.x_neurons x];
        Neurons.y_neurons = [Neurons.y_neurons y];
        Neurons.z_neurons = [Neurons.z_neurons z];
        for j = 1 : Electrode.nchannels
            distance_aux(j,1) = sqrt((Neurons.x_neurons(1,end)-Electrode.coordinates(j,1)).^2+(Neurons.y_neurons(1,end)-Electrode.coordinates(j,2)).^2+...
            (Neurons.z_neurons(1,end)-Electrode.coordinates(j,3)).^2);
        end
        distance = [distance distance_aux];
    end
end

% Classify between close and distant neurons
%distance
%Par_sim.dlim
%Electrode.diameter
%Electrode.nchannels
%disp('distance==');
%distance 
disp([' Par_sim.dlim=' num2str(Par_sim.dlim) ...
      ' Electrode.diameter=' num2str(Electrode.diameter) ' Electrode.nchannels=' num2str(Electrode.nchannels)]);
volume_close = ceil(find(distance<=(Par_sim.dlim+Electrode.diameter/2))./Electrode.nchannels);   % Neurons to be simulated in detail
%disp('volume_close=');
%volume_close
repeated = diff(volume_close)==0;
volume_close(repeated) = [];
volume_far = 1:length(Neurons.x_neurons);
%disp('volume_far  =');
%volume_far
volume_far(volume_close) = [];
nneurons_close = round(length(volume_close)); 
nneurons_far = length(volume_far);

nneurons_far_plot = round(nneurons_far*Par_cube.prop_plot/100);
if Par_sim.orig_code
    nneurons_close_plot = round(nneurons_close*Par_cube.prop_plot/100);
else
    nneurons_close_plot = nneurons_close;
end
disp(['nneurons_far_plot=' num2str(nneurons_far_plot) ' nneurons_close_plot=' num2str(nneurons_close_plot) ]);

Neurons.nneurons = nneurons_close+nneurons_far;
Neurons.coordinates = [];
Neurons.coordinates = [Neurons.x_neurons' Neurons.y_neurons' Neurons.z_neurons'];

if Par_sim.manual == 1
    volume_close_auto = ceil(find(distance_auto<=(Par_sim.dlim+Electrode.diameter/2))./Electrode.nchannels);   % Neurons to be simulated in detail
    repeated = diff(volume_close_auto)==0;
    volume_close_auto(repeated) = [];
    volume_far_auto = 1:length(x_neurons_auto);
    volume_far_auto(volume_close_auto) = [];
    nneurons_close_auto = round(length(volume_close_auto)); 
    nneurons_far_auto = length(volume_far_auto);
    Neurons_aux.nneurons = nneurons_close_auto + nneurons_far_auto;
    Neurons_aux.coordinates = [];
    Neurons_aux.coordinates = [x_neurons_auto' y_neurons_auto' z_neurons_auto'];
end

% Plot the cube

cla(handles.Cube_plot)

%if Par_sim.orig_code
    x_neurons_far = Neurons.x_neurons(volume_far(1:nneurons_far_plot));
    y_neurons_far = Neurons.y_neurons(volume_far(1:nneurons_far_plot));
    z_neurons_far = Neurons.z_neurons(volume_far(1:nneurons_far_plot));
    plot3(handles.Cube_plot,x_neurons_far,y_neurons_far,z_neurons_far,'^c','markersize',3,'linewidth',3)
    set(handles.Cube_plot,'XGrid','on','YGrid','on','ZGrid','on')
    hold(handles.Cube_plot,'on')
%end

x_neurons_close = Neurons.x_neurons(volume_close(1:nneurons_close_plot));
y_neurons_close = Neurons.y_neurons(volume_close(1:nneurons_close_plot));
z_neurons_close = Neurons.z_neurons(volume_close(1:nneurons_close_plot));
plot3(handles.Cube_plot,x_neurons_close,y_neurons_close,z_neurons_close,'^r','markersize',3,'linewidth',3)
set(handles.Cube_plot,'XGrid','on','YGrid','on','ZGrid','on')
axis(handles.Cube_plot,[-distance_max distance_max -distance_max distance_max -distance_max distance_max])
hold(handles.Cube_plot,'on')

% Plot the electrode

set(handles.neurocube_figure,'CurrentAxes',handles.Cube_plot)
for i = 1 : Electrode.nchannels
    h = filledCircle(Electrode.coordinates(i,:),Electrode.diameter/2,100,'k');
end

% Assign neuron models to close neurons

Models = randi(5,1,nneurons_close);     % 1 = 'd151'; 2 = 'pyramidal264L'; 3 = '51-2a'; 4 = '51-2b'; 5 = 're80'
intern = find(Models==5);               % Interneurons
Parameters = randi(4,1,nneurons_close);
Neurons.id =[];
Neurons.id(volume_close) = 600 + (Models-1)*4 + Parameters;

if Par_sim.manual == 1
    Models_auto = randi(5,1,nneurons_close_auto);     % 1 = 'd151'; 2 = 'pyramidal264L'; 3 = '51-2a'; 4 = '51-2b'; 5 = 're80'
    intern_auto = find(Models_auto==5);               % Interneurons
    Parameters_auto = randi(4,1,nneurons_close_auto);
    Neurons_aux.id =[];
    Neurons_aux.id(volume_close_auto) = 600 + (Models_auto-1)*4 + Parameters_auto;
end

% Assign neuron models to distant neurons

load spike_shapes.mat
SpikeShapes = NormShapes./repmat(max(NormShapes,[],2),1,size(NormShapes,2));
clear Shapes
[nshapes,~] = size(SpikeShapes);
Neurons.id(volume_far) = randi(nshapes,1,nneurons_far);
if Par_sim.manual == 1
    Neurons_aux.id(volume_far_auto) = randi(nshapes,1,nneurons_far_auto);
end

% Assign firing rates

if Par_sim.orig_firing_rates
    if Par_cube.firing_rates == 1   % Exponential
        k_exp = Par_cube.firing_rates_params(1);
        sigma_exp = Par_cube.firing_rates_params(2);
        theta_exp = Par_cube.firing_rates_params(3);
        Rates = gprnd(k_exp,sigma_exp,theta_exp,1,Neurons.nneurons);
        high_freq = find(Rates>10);
        Rates(high_freq) = 0.5*rand(1,length(high_freq));
        if Par_sim.manual == 1
            Rates_auto = gprnd(k_exp,sigma_exp,theta_exp,1,Neurons_aux.nneurons);
            high_freq_auto = find(Rates_auto>10);
            Rates_auto(high_freq_auto) = 0.5*rand(1,length(high_freq_auto));
        end
    elseif Par_cube.firing_rates == 2   % Uniform
        min_rate = Par_cube.firing_rates_params(4);
        max_rate = Par_cube.firing_rates_params(5);
        Rates = min_rate + (max_rate-min_rate).*rand(Neurons.nneurons,1);
        if Par_sim.manual ==1
            Rates_auto = min_rate + (max_rate-min_rate).*rand(Neurons_aux.nneurons,1);
        end
    end
    Rates(volume_close(intern)) = Rates(volume_close(intern))*5;    % We increase the firing rate of close interneurons to make it 5 times larger
    Far_intern_aux = randperm(nneurons_far) > (nneurons_far*0.8);
    Far_intern = find(Far_intern_aux == 1);
    Rates(Far_intern) = Rates(Far_intern)*5;                        % We increase the firing rate of distant interneurons to make it 5 times larger
    if Par_sim.manual == 1
        for i = 1 : length(Neurons.nmanual)
            single_aux = Neurons.nmanual(i);
            real_rate = Neurons.rates_manual(single_aux);
            Rates(end-length(Neurons.nmanual)+i) = real_rate;
        end
        Rates_auto(volume_close_auto(intern_auto)) = Rates_auto(volume_close_auto(intern_auto))*5;  
        Far_intern_aux_auto = randperm(nneurons_far_auto) > (nneurons_far_auto*0.8);
        Far_intern_auto = find(Far_intern_aux_auto == 1);
        Rates_auto(Far_intern_auto) = Rates_auto(Far_intern_auto)*5;
        Neurons_aux.Rates = [];
        Neurons_aux.Rates = Rates_auto;
    end
else
    Rates = 5.* ones(1, Neurons.nneurons); % mean # of spikes per second
end
disp(' ');
disp('Rates=');
Rates
disp(' ');

Neurons.Rates = [];
Neurons.Rates = Rates;
Neurons.firing_rates = Rates;

% Calculate spike times for close neurons

tres = 1/(Par_sim.downsampling*Par_sim.sr);
Refract = Par_sim.refract_period * 0.001 / tres;
Neurons.spiketimes = [];
total_spikes_close_neurons = 0;
for i = 1 : nneurons_close
    ISI = exprnd(Par_sim.sr*Par_sim.downsampling/Neurons.Rates(volume_close(i)),1,round(Neurons.Rates(volume_close(i))*Par_sim.duration));
    spiketimes = round(cumsum(ISI));
    CloseInds = diff(spiketimes) <= Refract; % Exclude potential overlapping spikes for a given neuron
    spiketimes(CloseInds) = [];
    LateInds = spiketimes>Par_sim.duration*Par_sim.sr*Par_sim.downsampling;
    spiketimes(LateInds) = [];
    Neurons.spiketimes{volume_close(i)} = spiketimes;
    disp(['Close neuron' num2str(i) ': num_spikes =' num2str(numel(spiketimes))]);
    total_spikes_close_neurons = total_spikes_close_neurons + numel(spiketimes);
end
disp(' ');
disp(['total_spikes_close_neurons=' num2str(total_spikes_close_neurons)]);

if Par_sim.manual ==1
    Neurons_aux.spiketimes =[];
    for i = 1 : nneurons_close_auto
        ISI_auto = exprnd(Par_sim.sr*Par_sim.downsampling/Neurons_aux.Rates(volume_close_auto(i)),1,round(Neurons_aux.Rates(volume_close_auto(i))*Par_sim.duration));
        spiketimes_auto = round(cumsum(ISI_auto));
        CloseInds_auto = diff(spiketimes_auto) <= Refract; % Exclude potential overlapping spikes for a given neuron
        spiketimes_auto(CloseInds_auto) = [];
        LateInds_auto = spiketimes_auto>Par_sim.duration*Par_sim.sr*Par_sim.downsampling;
        spiketimes_auto(LateInds_auto) = [];
        Neurons_aux.spiketimes{volume_close_auto(i)} = spiketimes_auto;
    end
end

USER_DATA{1} = Par_cube;  
USER_DATA{2} = Neurons;  
USER_DATA{3} = Electrode;  
USER_DATA{4} = Par_sim;
set(handles.neurocube_figure,'userdata',USER_DATA);

if Par_sim.manual == 1
    Neurons_aux.nmanual = [];
    Neurons_aux.distance_manual = Neurons.distance_manual;        
    Neurons_aux.rates_manual = Neurons.rates_manual;
    USER_DATA{5} = Neurons_aux;
else    % Auto mode
    Neurons_temp = Neurons;
    Par_sim.manual = 1;
    USER_DATA{4} = Par_sim;
    set(handles.neurocube_figure,'userdata',USER_DATA);
    Update_auto(handles);
    USER_DATA = get(handles.neurocube_figure,'userdata');
    USER_DATA{5} = USER_DATA{2};
    USER_DATA{2} = Neurons_temp;
    Par_sim = USER_DATA{4};
    Par_sim.manual = 0;
    USER_DATA{4} = Par_sim;
    set(handles.neurocube_figure,'userdata',USER_DATA);
end

set(handles.Progress_text,'Visible','Off');
set(handles.Save_cube,'Enable','on')
set(handles.Run,'Enable','on')
set(handles.Update_cube,'Enable','on')

set(handles.neurocube_figure,'userdata',USER_DATA);

function Cube_edge_Callback(hObject, eventdata, handles)
% hObject    handle to Cube_edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cube_edge as text
%        str2double(get(hObject,'String')) returns contents of Cube_edge as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Par_cube.edge = str2double(get(hObject,'String'));
USER_DATA{1} = Par_cube;
set(handles.neurocube_figure,'userdata',USER_DATA);

% --- Executes during object creation, after setting all properties.
function Cube_edge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cube_edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Par_cube.density = str2double(get(hObject,'String'));
USER_DATA{1} = Par_cube;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prop_firing_Callback(hObject, eventdata, handles)
% hObject    handle to Prop_firing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prop_firing as text
%        str2double(get(hObject,'String')) returns contents of Prop_firing as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Par_cube.prop_firing = str2double(get(hObject,'String'));
USER_DATA{1} = Par_cube;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function Prop_firing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prop_firing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in Firing_rates.
function Firing_rates_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Firing_rates 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
if eventdata.NewValue == handles.Exp
    Par_cube.firing_rates = 1;
    set(handles.Min_freq,'Enable','off')
    set(handles.Max_freq,'Enable','off')
elseif eventdata.NewValue == handles.Unif
    Par_cube.firing_rates = 2;
    set(handles.Min_freq,'Enable','on')
    set(handles.Max_freq,'Enable','on')
end
USER_DATA{1} = Par_cube;
set(handles.neurocube_figure,'userdata',USER_DATA);



function Duration_sim_Callback(hObject, eventdata, handles)
% hObject    handle to Duration_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Duration_sim as text
%        str2double(get(hObject,'String')) returns contents of Duration_sim as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_sim = USER_DATA{4};
Par_sim.duration = str2double(get(hObject,'String'));
USER_DATA{4} = Par_sim;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function Duration_sim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Duration_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load_cube.
function Load_cube_Callback(hObject, eventdata, handles)
% hObject    handle to Load_cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.neurocube_figure,'userdata');
[file,path] = uigetfile('*.mat','Load file name');
file_name = [path filesep file];
load(file_name)
USER_DATA{2} = Neurons;
Par_sim = USER_DATA{4};
if ~isempty(Neurons.nmanual)
    set(handles.Auto,'Value',0);
    set(handles.Manual,'Value',1);
    Par_sim.manual = 1;
else
    set(handles.Auto,'Value',1);
    set(handles.Manual,'Value',0);
    Par_sim.manual = 0;
end

if exist('Neurons_aux','var')
    USER_DATA{5} = Neurons_aux;
end

USER_DATA{4} = Par_sim;

set(handles.Save_cube,'Enable','on')
set(handles.Run,'Enable','on')
set(handles.Update_cube,'Enable','on')

set(handles.neurocube_figure,'userdata',USER_DATA);
Update_plot(handles)



% --- Executes on button press in Save_cube.
function Save_cube_Callback(hObject, eventdata, handles)
% hObject    handle to Save_cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
if length(USER_DATA) > 4
    Neurons_aux = USER_DATA{5};
    [file,path] = uiputfile('Neurons.mat','Save file name');
    file_name = [path filesep file];
    save(file_name,'Neurons','Neurons_aux')
else
    [file,path] = uiputfile('Neurons.mat','Save file name');
    file_name = [path filesep file];
    save(file_name,'Neurons')
end


% --- Executes when selected object is changed in Nchannels.
function Nchannels_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Nchannels 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.neurocube_figure,'userdata');
Electrode = USER_DATA{3};
Neurons = USER_DATA{2};
if eventdata.NewValue == handles.Single
    Electrode.nchannels = 1;
    Electrode.coordinates = [0 0 0];
    set(handles.Separation,'Enable','off')
elseif eventdata.NewValue == handles.Tetrode
    Electrode.nchannels = 4;
    Electrode.coordinates = (Electrode.separation/2)*Electrode.tetrode;
    set(handles.Separation,'Enable','on')
end
USER_DATA{3} = Electrode;
set(handles.neurocube_figure,'userdata',USER_DATA);
if ~isempty(Neurons.coordinates)
    Update_neurons(handles)
    Update_plot(handles)
end



function Diameter_Callback(hObject, eventdata, handles)
% hObject    handle to Diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Diameter as text
%        str2double(get(hObject,'String')) returns contents of Diameter as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Electrode = USER_DATA{3};
Neurons = USER_DATA{2};
Electrode.diameter = str2double(get(hObject,'String'));
USER_DATA{3} = Electrode;
set(handles.neurocube_figure,'userdata',USER_DATA);
if ~isempty(Neurons.coordinates)
    Update_neurons(handles)
    Update_plot(handles)
end


% --- Executes during object creation, after setting all properties.
function Diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Separation_Callback(hObject, eventdata, handles)
% hObject    handle to Separation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Separation as text
%        str2double(get(hObject,'String')) returns contents of Separation as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Electrode = USER_DATA{3};
Neurons = USER_DATA{2};
Separation = str2double(get(hObject,'String'));
if Electrode.nchannels == 1
    Electrode.coordinates = [0 0 0];
elseif Electrode.nchannels == 4
    Electrode.coordinates = (Separation/2)*Electrode.tetrode;
end
Electrode.separation = Separation;
USER_DATA{3} = Electrode;
set(handles.neurocube_figure,'userdata',USER_DATA);
if ~isempty(Neurons.coordinates)
    Update_neurons(handles)
    Update_plot(handles)
end


% --- Executes during object creation, after setting all properties.
function Separation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Separation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_simulation(handles);


% --- Executes on button press in Save_sim.
function Save_sim_Callback(hObject, eventdata, handles)
% hObject    handle to Save_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Par_sim = USER_DATA{4};
[file,path] = uiputfile('Spikesim.mat','Save file name');
file_name = [path filesep file];

set(handles.Progress_text,'String','Saving file...');
set(handles.Progress_text,'Visible','On');
pause(0.01)

data = Neurons.data;
Neurons_id = Neurons.id;
Coordinates = Neurons.coordinates;

volume_close = find(Neurons.id >= 600);

samples_per_ms = 1000/(Par_sim.sr*Par_sim.downsampling);
for i = 1 : length(Neurons.spiketimes)
    Spiketimes{i} = samples_per_ms.*Neurons.spiketimes{i};
end

x_neurons = Neurons.x_neurons;
y_neurons = Neurons.y_neurons;
z_neurons = Neurons.z_neurons;
models    = Neurons.models;
firing_rates = Neurons.firing_rates;
nneurons_far = Neurons.nneurons_far;

Close_Spikeshapes = Neurons.spikeshapes;
Coordinates_close = Coordinates(volume_close,:);
Spiketimes_close = Spiketimes(volume_close);
for i = 1 : length(Spiketimes_close)
    Spikeclass_close{i} = i*ones(size(Spiketimes_close{i}));
end
Close_neurons(:,1) = cell2mat(Spikeclass_close);
Close_neurons(:,2) = cell2mat(Spiketimes_close);
Close_neurons = sortrows(Close_neurons,2);


save(file_name,'data','Neurons_id','Spiketimes','Coordinates','Close_Spikeshapes',...
               'Coordinates_close','Close_neurons','x_neurons','y_neurons','z_neurons',...
               'models', 'firing_rates', 'nneurons_far', 'Par_sim')

set(handles.Progress_text,'String','File saved');



function Min_freq_Callback(hObject, eventdata, handles)
% hObject    handle to Min_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min_freq as text
%        str2double(get(hObject,'String')) returns contents of Min_freq as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Par_cube.firing_rates_params(4) = str2double(get(hObject,'String'));
USER_DATA{1} = Par_cube;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function Min_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Max_freq_Callback(hObject, eventdata, handles)
% hObject    handle to Max_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Max_freq as text
%        str2double(get(hObject,'String')) returns contents of Max_freq as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Par_cube.firing_rates_params(5) = str2double(get(hObject,'String'));
USER_DATA{1} = Par_cube;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function Max_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Progress_text_Callback(hObject, eventdata, handles)
% hObject    handle to Progress_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Progress_text as text
%        str2double(get(hObject,'String')) returns contents of Progress_text as a double


% --- Executes during object creation, after setting all properties.
function Progress_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Progress_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance1_Callback(hObject, eventdata, handles)
% hObject    handle to distance1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance1 as text
%        str2double(get(hObject,'String')) returns contents of distance1 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.distance_manual(1) = str2double(get(hObject,'String'));
if Neurons.distance_manual(1) < 0
    Neurons.distance_manual(1) = 0;
    set(handles.distance1,'String','0');
end
if Neurons.distance_manual(1) > 1
    Neurons.distance_manual(1) = 1;
    set(handles.distance1,'String','1');
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function distance1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function distance2_Callback(hObject, eventdata, handles)
% hObject    handle to distance2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance2 as text
%        str2double(get(hObject,'String')) returns contents of distance2 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.distance_manual(2) = str2double(get(hObject,'String'));
if Neurons.distance_manual(2) < 0
    Neurons.distance_manual(2) = 0;
    set(handles.distance2,'String','0');
end
if Neurons.distance_manual(2) > 1
    Neurons.distance_manual(2) = 1;
    set(handles.distance2,'String','1');
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function distance2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance3_Callback(hObject, eventdata, handles)
% hObject    handle to distance3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance3 as text
%        str2double(get(hObject,'String')) returns contents of distance3 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.distance_manual(3) = str2double(get(hObject,'String'));
if Neurons.distance_manual(3) < 0
    Neurons.distance_manual(3) = 0;
    set(handles.distance3,'String','0');
end
if Neurons.distance_manual(3) > 1
    Neurons.distance_manual(3) = 1;
    set(handles.distance3,'String','1');
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function distance3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance4_Callback(hObject, eventdata, handles)
% hObject    handle to distance4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance4 as text
%        str2double(get(hObject,'String')) returns contents of distance4 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.distance_manual(4) = str2double(get(hObject,'String'));
if Neurons.distance_manual(4) < 0
    Neurons.distance_manual(4) = 0;
    set(handles.distance4,'String','0');
end
if Neurons.distance_manual(4) > 1
    Neurons.distance_manual(4) = 1;
    set(handles.distance4,'String','1');
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function distance4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance5_Callback(hObject, eventdata, handles)
% hObject    handle to distance5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance5 as text
%        str2double(get(hObject,'String')) returns contents of distance5 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.distance_manual(5) = str2double(get(hObject,'String'));
if Neurons.distance_manual(5) < 0
    Neurons.distance_manual(5) = 0;
    set(handles.distance5,'String','0');
end
if Neurons.distance_manual(5) > 1
    Neurons.distance_manual(5) = 1;
    set(handles.distance5,'String','1');
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function distance5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rate1_Callback(hObject, eventdata, handles)
% hObject    handle to rate1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rate1 as text
%        str2double(get(hObject,'String')) returns contents of rate1 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.rates_manual(1) = str2double(get(hObject,'String'));
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function rate1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rate2_Callback(hObject, eventdata, handles)
% hObject    handle to rate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rate2 as text
%        str2double(get(hObject,'String')) returns contents of rate2 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.rates_manual(2) = str2double(get(hObject,'String'));
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function rate2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rate3_Callback(hObject, eventdata, handles)
% hObject    handle to rate3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rate3 as text
%        str2double(get(hObject,'String')) returns contents of rate3 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.rates_manual(3) = str2double(get(hObject,'String'));
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function rate3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rate4_Callback(hObject, eventdata, handles)
% hObject    handle to rate4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rate4 as text
%        str2double(get(hObject,'String')) returns contents of rate4 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.rates_manual(4) = str2double(get(hObject,'String'));
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function rate4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rate5_Callback(hObject, eventdata, handles)
% hObject    handle to rate5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rate5 as text
%        str2double(get(hObject,'String')) returns contents of rate5 as a double
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
Neurons.rates_manual(5) = str2double(get(hObject,'String'));
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes during object creation, after setting all properties.
function rate5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in neuron1.
function neuron1_Callback(hObject, eventdata, handles)
% hObject    handle to neuron1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuron1
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
selected = get(hObject,'Value');
if selected == 1
    Neurons.nmanual = sort([Neurons.nmanual 1]);
else
    i = Neurons.nmanual == 1;
    Neurons.nmanual(i) = [];
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes on button press in neuron2.
function neuron2_Callback(hObject, eventdata, handles)
% hObject    handle to neuron2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuron2
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
selected = get(hObject,'Value');
if selected == 1
    Neurons.nmanual = sort([Neurons.nmanual 2]);
else
    i = Neurons.nmanual == 2;
    Neurons.nmanual(i) = [];
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);

% --- Executes on button press in neuron3.
function neuron3_Callback(hObject, eventdata, handles)
% hObject    handle to neuron3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuron3
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
selected = get(hObject,'Value');
if selected == 1
    Neurons.nmanual = sort([Neurons.nmanual 3]);
else
    i = Neurons.nmanual == 3;
    Neurons.nmanual(i) = [];
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes on button press in neuron4.
function neuron4_Callback(hObject, eventdata, handles)
% hObject    handle to neuron4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuron4
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
selected = get(hObject,'Value');
if selected == 1
    Neurons.nmanual = sort([Neurons.nmanual 4]);
else
    i = Neurons.nmanual == 4;
    Neurons.nmanual(i) = [];
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes on button press in neuron5.
function neuron5_Callback(hObject, eventdata, handles)
% hObject    handle to neuron5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of neuron5
USER_DATA = get(handles.neurocube_figure,'userdata');
Neurons = USER_DATA{2};
selected = get(hObject,'Value');
if selected == 1
    Neurons.nmanual = sort([Neurons.nmanual 5]);
else
    i = Neurons.nmanual == 5;
    Neurons.nmanual(i) = [];
end
USER_DATA{2} = Neurons;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes when selected object is changed in Single_units.
function Single_units_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Single_units 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
USER_DATA = get(handles.neurocube_figure,'userdata');
Par_sim = USER_DATA{4};
if eventdata.NewValue == handles.Auto
    Par_sim.manual = 0;
    Neurons_aux = USER_DATA{2};
    Neurons = USER_DATA{5};
    USER_DATA{5} = Neurons_aux;
    USER_DATA{2} = Neurons;
    set(handles.distance1,'Enable','off')
    set(handles.distance2,'Enable','off')
    set(handles.distance3,'Enable','off')
    set(handles.distance4,'Enable','off')
    set(handles.distance5,'Enable','off')
    set(handles.rate1,'Enable','off')
    set(handles.rate2,'Enable','off')
    set(handles.rate3,'Enable','off')
    set(handles.rate4,'Enable','off')
    set(handles.rate5,'Enable','off')
    set(handles.neuron1,'Enable','off')
    set(handles.neuron2,'Enable','off')
    set(handles.neuron3,'Enable','off')
    set(handles.neuron4,'Enable','off')
    set(handles.neuron5,'Enable','off')
elseif eventdata.NewValue == handles.Manual
    Par_sim.manual = 1;
    Neurons_aux = USER_DATA{2};
    if length(USER_DATA) > 4
        Update_auto(handles)
        USER_DATA{2} = USER_DATA{5};
        Update_plot(handles)
    end
    USER_DATA{5} = Neurons_aux;
    set(handles.distance1,'Enable','on')
    set(handles.distance2,'Enable','on')
    set(handles.distance3,'Enable','on')
    set(handles.distance4,'Enable','on')
    set(handles.distance5,'Enable','on')
    set(handles.rate1,'Enable','on')
    set(handles.rate2,'Enable','on')
    set(handles.rate3,'Enable','on')
    set(handles.rate4,'Enable','on')
    set(handles.rate5,'Enable','on')
    set(handles.neuron1,'Enable','on')
    set(handles.neuron2,'Enable','on')
    set(handles.neuron3,'Enable','on')
    set(handles.neuron4,'Enable','on')
    set(handles.neuron5,'Enable','on')
end
USER_DATA{4} = Par_sim;
set(handles.neurocube_figure,'userdata',USER_DATA);


% --- Executes on button press in Update_cube.
function Update_cube_Callback(hObject, eventdata, handles)
% hObject    handle to Update_cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
Update_plot(handles)
