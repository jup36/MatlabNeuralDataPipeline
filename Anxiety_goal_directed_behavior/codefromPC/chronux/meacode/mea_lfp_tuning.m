function lfp_tuning=mea_lfp_tuning(no_electrode,spec_data,avg_flag)
%% lfp_tuning=mea_lfp_tuning(no_electrode)
if ~exist('no_electrode','var')
    disp('which electrode do you want to analyze?   ');
end

if ~exist('spec_data','var')
    try 
        spec_data=evalin('base','spec_data');
        disp('no input spectrum data. by default, use the data in base workspace');
    catch
        disp('no input spectrum data!');
        return;
    end
end

if ~exist('avg_flag','var')
    avg_flag=0;
end

[num_trial,num_direction]=size(spec_data);
temp=spec_data{1,1};
[num_freq,num_ch]=size(temp.freq_power);
if no_electrode>num_ch
    disp('no selected channel!');
    return;
end

spec_mat=zeros(num_freq,num_direction,num_trial);
for m1=1:num_trial
    for m2=1:num_direction
        temp=spec_data{m1,m2};
        spec_mat(:,m2,m1)=temp.freq_power(:,no_electrode);
    end
end

if avg_flag==1
    lfp_tuning=mean(spec_mat,3);
else
    lfp_tuning=spec_mat;
end

return;