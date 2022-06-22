function norm_neuron=normalize_spike_shape2(neuron,nsamples,nrow)
% Addapt the number of samples of the spike shape 'neuron' to nsamples

shape=neuron(nrow,4:end);
t=neuron(1,4:end);
step=(t(end)-t(1))/(nsamples-1);
t_interp=t(1):step:t(end);
shape_interp=interp1(t,shape,t_interp);
norm_shape=shape_interp;
norm_neuron(1,:)=t_interp;
norm_neuron(2,:)=norm_shape;