function [fire_rate,direction,sum_spike]=mea_tuning(events,no_neuron,params)

if ~exist('events','var')
    display('please select dataset!');
end
%in case no enough input variables
try 
    pd=evalin('base','pd');
    events=pd.EVENTS;
    clear pd;
catch
    display('no input dataset!');
    return;
end

if ~exist('no_neuron','var')||(no_neuron==0)
    display('you didn''t select neuron NOs! By default, all neurons will be analyzed'); 
    temp=size(events);
    no_neuron=1:temp(1);
end
%in case no neuron number assigned

try
    fs=params.Fs;
catch 
    fs=1000;
end
%sampling frequency

try 
    t_step=params.t_step;
catch
    t_step=input('period for each stimulus:     ');
end
%period for each stimulus

try 
    drift_t=params.drift_t;
catch
    drift_t=0;
end
%drift t

try 
    num_trial=params.num_trial;
catch %#ok<*CTCH>
    [~,num_trial]=size(events);
end
%get the number of trials

try 
    use_trial=params.use_trial;
catch
    use_trial=1:num_trial;
end
%use_trial choose trials that we need to calculate. 

try 
    direction=params.direction;
catch
    num_direction=input('how many directions do you use?:   ');
    direction=(1:num_direction)*2*pi/360;
end
% orientation stimulus 

num_spike=zeros(length(no_neuron),length(direction));
%used to store number of spikes for each direction

for m1=1:length(no_neuron)
    no1=no_neuron(m1);
    for m2=use_trial
        tmp_ts=events{no1,m2}+drift_t;
        if isempty(tmp_ts)
            continue;
        end
        %when time stamps are empty, jump out of the loop.
        temp=floor(uint32(tmp_ts*fs)/t_step);
        
        temp(temp>length(direction))=direction(end);
        [tmp_num,~]=hist(temp,1:length(direction));
        num_spike(m1,:)=num_spike(m1,:)+tmp_num;
    end
end
%calculate total number of spikes for each stimulus.

fire_rate=num_spike/(t_step/fs)/length(use_trial);
%unit of fire_rate: Hz
sum_spike=sum(num_spike,2);
%total number of spikes for each neuron

return;














