function run_simulation(handles)
% Calculate the extracellular recording generated by the neuron
% distribution

USER_DATA = get(handles.neurocube_figure,'userdata');
Par_cube = USER_DATA{1};
Neurons = USER_DATA{2};
Electrode = USER_DATA{3};
Par_sim = USER_DATA{4};

set(handles.Run,'Enable','off')
set(handles.Save_sim,'Enable','off')

x_neurons = Neurons.coordinates(:,1);
y_neurons = Neurons.coordinates(:,2);
z_neurons = Neurons.coordinates(:,3);
volume_close = find(Neurons.id >= 600);
volume_far = find(Neurons.id < 600);
Models = ceil((Neurons.id(volume_close)-600)./4);
if Par_sim.five_parameters
    Parameters = mod((Neurons.id(volume_close)-601),4)+1;
else
    num_neurons = numel(x_neurons);
    Parameters = 3 * ones(1, num_neurons);    
end    

if 1
    if Electrode.diameter < 7
        Electrode.diameter = 7;
        set(handles.Diameter,'String','7');
    end
end
rad = Electrode.diameter/2;
nchannels = Electrode.nchannels;
nneurons = Neurons.nneurons;

nsim_points = Par_sim.grid_npoint;
gridsize = round(Electrode.diameter/(sqrt(nsim_points)-1));

neuron_models = Par_sim.neuron_models;

disp(['Models=' cellstr(Par_sim.neuron_models') ' Parameters=' num2str(Parameters) ' volume_close=' num2str(volume_close) ' volume_far=' num2str(volume_far)]);

model_center = [0 0 0;-5 -5 0;-5 -10 -45;0 20 50;-10 0 20];

allshapes = zeros(nneurons*nchannels,316);

allamps = ones(1,nneurons*nchannels);

Shapes = zeros(nneurons*nchannels,316);

% Calculate the spike shapes for close neurons

set(handles.Progress_text,'String','Calculating spike shapes for close-by neurons...');
set(handles.Progress_text,'Visible','On');
pause(0.01)

nstep_bar = round(length(volume_close)/10);

if Par_sim.orig_code
    for i = 1 : length(volume_close)
        for j = 1 : nchannels
            neuron_position=[x_neurons(volume_close(i)) y_neurons(volume_close(i)) z_neurons(volume_close(i))];
            soma_coord=neuron_position+model_center(Models(i),:);
            electrode_rel=neuron_position-Electrode.coordinates(j,:);
            electrode_virtual=soma_coord+Rotate_coordinates(electrode_rel,Models(i));
            final_coord=electrode_virtual-neuron_position;
            x=final_coord(1);
            y=final_coord(2);
            z=final_coord(3);
            if rad == 0
                neuron = zplane_eap_calc(neuron_models{Models(i)}, Parameters(i), [], 0.3, [x x y y], z, gridsize);
                neuron_norm = 1000.*normalize_spike_shape(neuron,316);                  % Convert to microvolts
                allshapes((volume_close(i)-1)*nchannels+j,:) = real(neuron_norm(2,:));  
            else
                neuron = zplane_eap_calc(neuron_models{Models(i)}, Parameters(i), [], 0.3, [x-rad x+rad y-rad y+rad], z, gridsize);
                disp('i, j, final_coord=');
                [i j final_coord]
                disp(['size(neuron)=' num2str(size(neuron))]);
                
                for k = 2:size(neuron,1)
                    neuron_norm = 1000.*normalize_spike_shape2(neuron,316,k);
                    spike_aux(k-1,:) = real(neuron_norm(2,:));      
                    disp(['k=' num2str(k) ' size(neuron_norm)=' num2str(size(neuron_norm)) ' size(spike_aux)=' num2str(size(spike_aux))]);                    
                end
                max_aux = max(abs(spike_aux)');                                         % Eliminates NaN
                index_aux = max_aux<1e10;
                average_spike = mean(spike_aux(index_aux,:),1);
                disp([' size(average_spike)=' num2str(size(average_spike))]);
                allshapes((volume_close(i)-1)*nchannels+j,:) = average_spike;
                disp(['i=' num2str(i) ' j=' num2str(j) 'size(spike_aux)=' num2str(size(spike_aux)) ' size(average_spike)=' num2str(size(average_spike))]);    
                disp(['i=' num2str(i) ' j=' num2str(j) ' size(average_spike)=' num2str(size(average_spike))]);
            end
            if Models(i) == 4 %51-2b
                allshapes((volume_close(i)-1)*nchannels+j,:) = allshapes((volume_close(i)-1)*nchannels+j,:)./2;
            end
        end
        if mod(i,nstep_bar) == 0
            Update_progress(handles,i,length(volume_close))
        end
    end
else % Par_sim.orig_code == 0
    for i = 1 : length(volume_close)
        for j = 1 : nchannels
            neuron_position=[x_neurons(volume_close(i)) y_neurons(volume_close(i)) z_neurons(volume_close(i))];
            soma_coord=neuron_position+model_center(Models(i),:);
            neuron_position
            Electrode.coordinates(j,:)
            electrode_rel=neuron_position-Electrode.coordinates(j,:);
            electrode_virtual=soma_coord+Rotate_coordinates(electrode_rel,Models(i));
            final_coord=electrode_virtual-neuron_position;
            x=final_coord(1);
            y=final_coord(2);
            z=final_coord(3);
            disp('i, j, final_coord=');
            [i j final_coord]
            disp(' ');
            electrode_rel
%           my_y= electrode_rel * electrode_rel'
            dist = sqrt(electrode_rel * electrode_rel');
              
            disp(['i=' num2str(i) ' j=' num2str(j) ' dist=' num2str(dist) ' rel_position=' num2str(electrode_rel) ' neuron_position=' num2str(neuron_position) ' Electrode.coordinates=' num2str(Electrode.coordinates(j,:))]);
            dist = dist*1.e-6;
            disp(' ');
            if rad == 0
                neuron = calculate_ampl_spike1(neuron_models{Models(i)}, Parameters(i), [], 0.3, [x x y y], z, gridsize, dist, Par_sim);
                neuron_norm = 1000.*normalize_spike_shape(neuron,316);                  % Convert to microvolts
                allshapes((volume_close(i)-1)*nchannels+j,:) = real(neuron_norm(2,:));
            else
                neuron = calculate_ampl_spike1(neuron_models{Models(i)}, Parameters(i), [], 0.3, [x-rad x+rad y-rad y+rad], z, gridsize, dist, Par_sim);
                for k = 2:size(neuron,1)
                    neuron_norm = 1000.*normalize_spike_shape2(neuron,316,k);
                    spike_aux(k-1,:) = real(neuron_norm(2,:));
                end
                max_aux = max(abs(spike_aux)');                                         % Eliminates NaN
                index_aux = max_aux<1e10;
                average_spike = mean(spike_aux(index_aux,:),1);
                allshapes((volume_close(i)-1)*nchannels+j,:) = average_spike;
                disp(['i=' num2str(i) ' j=' num2str(j) 'size(spike_aux)=' num2str(size(spike_aux)) ' size(average_spike)=' num2str(size(average_spike))]);
                disp(['i=' num2str(i) ' j=' num2str(j) ' size(average_spike)=' num2str(size(average_spike))]);
            end
            if Models(i) == 4 %51-2b
                allshapes((volume_close(i)-1)*nchannels+j,:) = allshapes((volume_close(i)-1)*nchannels+j,:)./2;
            end
        end
        if mod(i,nstep_bar) == 0
            Update_progress(handles,i,length(volume_close))
        end
    end
end

Neurons.spikeshapes = [];
for i = 1 : nchannels
    Neurons.spikeshapes(:,(i-1)*316+1:(i-1)*316+316) = allshapes((volume_close(:)-1)*nchannels+i,:);
end

Neurons.models = {};
for i = 1 : length(volume_close)
    Neurons.models{i} = neuron_models{Models(i)};
end

set(handles.Progress_text,'Visible','Off');
set(handles.Progress,'Visible','Off');
cla(handles.Progress)

% Calculate the spike shapes for distant neurons

set(handles.Progress_text,'String','Calculating spike shapes for distant neurons...');
set(handles.Progress_text,'Visible','On');
pause(0.01)

load spike_shapes.mat
SpikeShapes = NormShapes./repmat(max(NormShapes,[],2),1,size(NormShapes,2));
clear Shapes

load Model_amplitudes.mat
[~,index]=min(abs(positions-Par_sim.dlim));
k = mean_amplitude(index)/model2(index-1);

grid_i=-rad:gridsize:rad;
grid_j=-rad:gridsize:rad;
for i=1:length(grid_i)
shift_matrix(1,(i-1)*length(grid_j)+1:(i-1)*length(grid_j)+length(grid_i))=grid_i(i)*ones(1,length(grid_j));
shift_matrix(2,(i-1)*length(grid_j)+1:(i-1)*length(grid_j)+length(grid_i))=grid_j(:);
end
shift_matrix(3,:)=0;

nstep_bar = round(length(volume_far)/10);

for i = 1 : length(volume_far)
    shape_aux = -SpikeShapes(Neurons.id(volume_far(i)),:);
    for j = 1 : nchannels
%         x = -x_neurons(volume_far(i))+Electrode.coordinates(j,1);
%         y = -y_neurons(volume_far(i))+Electrode.coordinates(j,2);
%         z = -z_neurons(volume_far(i))+Electrode.coordinates(j,3);
        neuron_coord = [x_neurons(volume_far(i)) y_neurons(volume_far(i)) z_neurons(volume_far(i))];
        electrode_coord = [Electrode.coordinates(j,1) Electrode.coordinates(j,2) Electrode.coordinates(j,3)];
        allamps((volume_far(i)-1)*nchannels+j) = calculate_ampl_spike2(k,rad,neuron_coord,electrode_coord,shift_matrix);
        allshapes((volume_far(i)-1)*nchannels+j,:) = shape_aux;
    end
    
    if mod(i,nstep_bar) == 0
        Update_progress(handles,i,length(volume_far))
    end
end

set(handles.Progress_text,'Visible','Off');
set(handles.Progress,'Visible','Off');
cla(handles.Progress)

Shapes=bsxfun(@times,allamps',allshapes);

spikelen = size(Shapes,2);

% Calculate the spike times for distant neurons

set(handles.Progress_text,'String','Calculating spike times...');
set(handles.Progress_text,'Visible','On');
pause(0.01)

tres = 1/(Par_sim.downsampling*Par_sim.sr);
Ns1 = round(Par_sim.duration / tres);
Refract = Par_sim.refract_period * 0.001 / tres;

for i = 1 : length(volume_far)
    ISI = exprnd(Par_sim.sr*Par_sim.downsampling/Neurons.Rates(volume_far(i)),1,round(Neurons.Rates(volume_far(i))*Par_sim.duration));
    spiketimes = round(cumsum(ISI));
    CloseInds = diff(spiketimes) <= Refract;        % Exclude potential overlapping spikes for a given neuron
    spiketimes(CloseInds) = [];
    LateInds = spiketimes>Par_sim.duration*Par_sim.sr*Par_sim.downsampling;
    spiketimes(LateInds) = [];
    Neurons.spiketimes{volume_far(i)} = spiketimes;
    
    if mod(i,nstep_bar) == 0
        Update_progress(handles,i,length(volume_far))
    end
end

set(handles.Progress_text,'Visible','Off');
set(handles.Progress,'Visible','Off');
cla(handles.Progress)

% Segment data for long simulations 

max_segment_time_sec = 10; 
full_segments_number = floor(Par_sim.duration / max_segment_time_sec); % the total time is filled with this simulation segments. CPG
non_full_segment_time_sec = Par_sim.duration - (full_segments_number * max_segment_time_sec);
if full_segments_number > 0
    Ns1_i = round(max_segment_time_sec / tres);
else
    Ns1_i = 0;
end;
Ns1_f = Ns1 - (Ns1_i*full_segments_number);

% Create spiketimes, spikeclass and spikeamps

Neurons.spikeclass = Neurons.spiketimes;
for i = 1 : length(Neurons.spikeclass)
    Neurons.spikeclass{i}(:) = i;
end
spiketimes = cell2mat(Neurons.spiketimes);
spikeclass = cell2mat(Neurons.spikeclass);
spikeamps = ones(length(spiketimes),1);
[~, sortind] = sort(spiketimes);
spiketimes = spiketimes(sortind);
spikeclass = spikeclass(sortind);
spikeamps = spikeamps(sortind); 

set(handles.Progress_text,'String','Bulding final recording...');
set(handles.Progress_text,'Visible','On');
pause(0.01)

data = [];

for segment_i = 1 : full_segments_number
    % select the spikes corresponding to this segment
    [segment_spiketimes_i] = find (((segment_i-1)*max_segment_time_sec/tres)< spiketimes & ...
        spiketimes <= (segment_i*max_segment_time_sec/tres));
    spiketimes_i = spiketimes(segment_spiketimes_i)- ((segment_i-1)*max_segment_time_sec/tres); % spike times for this segment
    spikeclass_i = spikeclass(segment_spiketimes_i); % Spike class of the spikes in this segment
    spikeamps_i = spikeamps(segment_spiketimes_i); % Spike amplitudes of the spikes in this segment
    
    % Make the data vectors:

    [data_i] = zeros(Ns1_i,1);
    data_i_aux = data_i';
    clear data_i
    for i = 1 : nchannels
        data_i(i,:) = data_i_aux;
    end
    
    % Add spikes to vector
    for Si = 1 : length(spiketimes_i)
        this_samp  = spiketimes_i(Si);
        if (this_samp + spikelen) < Ns1_i
            for i= 1 : nchannels
                this_shape = spikeamps_i(Si).*Shapes((spikeclass_i(Si)-1)*nchannels+i,:);
                data_i(i,this_samp:this_samp+spikelen-1) = data_i(i,this_samp:this_samp+spikelen-1) + this_shape; % add spike to noise
            end
        end
     end
    %Downsample the data
    data_i = data_i(:,1:Par_sim.downsampling:length(data_i));
    data = [data data_i];
    clear data_i
end;

% Last segment - eventually shorter than the other ones
if non_full_segment_time_sec > 0
    % select the spikes corresponding to this segment
    [segment_spiketimes_i] = find ((full_segments_number*max_segment_time_sec/tres)< spiketimes & ...
        spiketimes <= ((full_segments_number*max_segment_time_sec+non_full_segment_time_sec)/tres));
    spiketimes_i = spiketimes(segment_spiketimes_i) - (full_segments_number*max_segment_time_sec/tres); % spike times for this segment
    spikeclass_i = spikeclass(segment_spiketimes_i); % Spike class of the spikes in this segment
    spikeamps_i = spikeamps(segment_spiketimes_i); % Spike amplitudes of the spikes in this segment
    
    % Make the data vectors:

    [data_i] = zeros(Ns1_f,1);
    data_i_aux = data_i';
    clear data_i
    for i = 1:nchannels
        data_i(i,:) = data_i_aux;
    end
    
    % Add spikes to vector
    for Si = 1 : length(spiketimes_i)
        this_samp  = spiketimes_i(Si);
        if (this_samp + spikelen) < Ns1_f
            for i = 1 : nchannels
                this_shape = spikeamps_i(Si).*Shapes((spikeclass_i(Si)-1)*nchannels+i,:);
                data_i(i,this_samp:this_samp+spikelen-1) = data_i(i,this_samp:this_samp+spikelen-1) + this_shape;
            end
        end
    end
    % downsample the data
    data_i = data_i(:,1:Par_sim.downsampling:length(data_i));
    data = [data data_i]; 
    clear data_i
end

set(handles.Progress_text,'Visible','Off');

% Filter the data

fmin = 300;
fmax = 3000;
sr = Par_sim.sr;
if exist('ellip','file') 
    [b,a] = ellip(2,0.1,40,[fmin fmax]*2/(sr));
else                % Filter coefficients for fmin=300, fmax=3000 and sr=24000
    a = [1.0000 -2.3930 2.0859 -0.9413 0.2502];
    b = [0.1966 -0.0167 -0.3598 -0.0167 0.1966];
end
length_data=size(data,2);
if Par_sim.orig_code
    th_noise=Par_sim.std_th_noise.*randn(1,length_data); 
end
[a_elec,b_elec] = filt_elec(2*rad);

% added by GD:
max_amplitude = [];
mean_amplitude = [];
max_amplitude_filtered = [];
mean_amplitude_filtered = [];
for i=1:nchannels
    if ~Par_sim.orig_code
        th_noise=Par_sim.std_th_noise.*randn(1,length_data);
    end
    max_amplitude = [max_amplitude  max(abs(data(i,:)))];
    mean_amplitude = [mean_amplitude mean(abs(data(i,:)))];
    if Par_sim.filter_data
        if exist('filtfilt','file')
            data_filt_elec(i,:) = filtfilt(b_elec,a_elec,data(i,:));
            data_noise(i,:) = data_filt_elec(i,:)+th_noise;
            data_filteredR(i,:) = filtfilt(b,a,data_noise(i,:));    % Final data with Robinson filter
        else
            data_filt_elec(i,:) = filter(b_elec,a_elec,data(i,:));
            data_noise(i,:) = data_filt_elec(i,:)+th_noise;
            data_filteredR(i,:) = filter(b,a,data_noise(i,:));    % Final data with Robinson filter
        end
        max_amplitude_filtered = [max_amplitude_filtered  max(abs(data_filteredR(i,:)))];
        mean_amplitude_filtered = [mean_amplitude_filtered  mean(abs(data_filteredR(i,:)))];
    end
end
data_initial=data;
if Par_sim.filter_data
    data=data_filteredR;
end
Neurons.data = data;

time_plot = 30;
t = 1/Par_sim.sr:1/Par_sim.sr:time_plot;
data_plot = data(:,1:length(t));

% Plot final data

cla(handles.Data_plot)

if Electrode.nchannels == 1
    set(handles.Data_plot,'Visible','On');
    plot(handles.Data_plot,t,data_plot)
else
    set(handles.Data_plot,'Visible','On');
    margin = 5;
    shift1 = max(data(1,:))+abs(min(data(2,:)))+margin;
    shift2 = max(data(2,:))+abs(min(data(3,:)))+margin;
    shift3 = max(data(3,:))+abs(min(data(4,:)))+margin;
    plot(handles.Data_plot,t,data_plot(1,:))
    hold(handles.Data_plot,'on')
    plot(handles.Data_plot,t,data_plot(2,:)+shift1)
    hold(handles.Data_plot,'on')
    plot(handles.Data_plot,t,data_plot(3,:)+shift1+shift2)
    hold(handles.Data_plot,'on')
    plot(handles.Data_plot,t,data_plot(4,:)+shift1+shift2+shift3)
end

disp(' ');
disp(['max_amplitude=' num2str(max_amplitude) ' length(volume_close)=' num2str(length(volume_close))]);
disp(['mean_amplitude=' num2str(mean_amplitude) ]);
if Par_sim.filter_data
    disp(['max_amplitude_filtered=' num2str(max_amplitude_filtered) ]);
    disp(['mean_amplitude_filtered=' num2str(mean_amplitude_filtered) ]);
end

for i=1:length(volume_close)
    num_spikes = numel(Neurons.spiketimes{volume_close(i)});
    disp(['cell#' num2str(i) ': num_spikes=' num2str(num_spikes)]);
end
disp(' ');
USER_DATA{2} = Neurons;
USER_DATA{3} = Electrode;

set(handles.Save_sim,'Enable','on')
set(handles.Run,'Enable','on')
set(handles.Progress_text,'String','Simulation ready');
set(handles.Progress_text,'Visible','On');

set(handles.neurocube_figure,'userdata',USER_DATA);
