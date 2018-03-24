function [pd pd_info]=mea_load(data_nm)
if nargin
    temp=load(data_nm);
else
    [tmp_file,tmp_dir]=uigetfile('*.mat');
    temp=load([tmp_dir,tmp_file]);
end;

tmp_nm=fieldnames(temp);
pd=getfield(temp,tmp_nm{1}); %#ok<GFLD>

[num_neuron,num_trial]=size(pd.EVENTS);
pd_info.num_neuron=num_neuron;%number of neurons
pd_info.num_trial=num_trial;%number of trials. 

[num_channel,total_t]=size(pd.NS2{1});
pd_info.num_channel=num_channel; %total number of channels
pd_info.total_t=total_t; %duration of time

pd_info.use_trial=1:num_trial;
return;
