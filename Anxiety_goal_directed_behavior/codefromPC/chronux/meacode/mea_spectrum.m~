function [spec_data,spec_info]=mea_spectrum(lfp,no_electrode,par)
%% [spec_data,spec_info]=mea_spectrum(lfp,no_electrode,par)
if ~exist('lfp','var')
    display('no input dataset, by default you will analyze variable pd in base workspace');
    try
        pd=evalin('base','pd');
        lfp=pd.NS2;
        clear pd;
        num_trial=length(lfp);
        [num_electrode total_time]=size(lfp{1});
        no_electrode=1:num_electrode;
    catch
        display('no enough input, quit!');
        return;
    end
else
    num_trial=length(lfp);
    [num_electrode total_time]=size(lfp{1});
    
end
% get dataset


if ~exist('no_electrode','var')
    display('you didn''t assign electrode to be analyzed, by defalt you will analyze all electrode');
    no_electrode=1:num_electrode;
elseif no_electrode==0
    no_electrode=1:num_electrode;
end

if isempty(no_electrode)
    display('no electrode selected');
    return;
end
% get electrode number to be analyzed.

try
    params.Fs=par.Fs;
catch
    params.Fs=input('sampling frequency:     ');
end
%sampling frequencyarams.

try
    params.fpass=par.fpass;
catch
    warning('missing frequency band pass is set to default [0 100]');
    params.fpass=[0 100];
end
% frequency band

try
    params.tapers=par.tapers;
catch
    warning('missing tapers is set to default [3 5]');
    params.tapers=[3 5];
end
% multitapers.

try
    params.err=par.err;
    err_flag=1;
catch
    temp=input('do you want to get confidence interval?(y/n)    :','s');
    if temp=='y'
        temp=input('error interval:     ','s');
        params.err=str2num(temp);
        err_flag=1;
    else
        err_flag=0;
    end
end

try
    t_start=par.t_start;
catch
    temp=input('time start:  ','s');
    t_start=str2num(temp);
end
% start time point of interval

try
    t_step=par.t_step;
catch
    temp=input('period for each interval:  ','s');
    t_step=str2num(temp);
end
t_end=t_start+t_step;
ind=t_end>total_time;
t_start(ind)=[];
t_end(ind)=[];
% time ending

if isempty(t_start)
    display('no interval to be analyzed');
    return;
end

try
    use_trial=params.use_trial;
catch
    display('you didn''t assign trials to be analyzed, by default it analyze all trials!');
    use_trial=1:num_trial;
end
% trials to be analyzed.

try
    mwin=params.mwin;
catch
    display('you didn''t assign value of moving window, by default is is set to [0.1 0.005]');
    mwin=[0.1,0.05];
end

if max(use_trial)>length(lfp)
    use_trial=1:length(lfp);
end

num_usetrial=length(use_trial);
num_useelec=length(no_electrode);
num_itval=length(t_start);
spec_data=cell(num_usetrial,num_itval);

for m1=1:num_usetrial
    no_trial=use_trial(m1);
    tmp_lfp=lfp{no_trial};
    temp=locdetrend(double(tmp_lfp'),params.Fs,mwin);
    tmp_lfp=temp';
    clear temp;
    if mod(uint16(m1*100/num_usetrial),10)==0
        display(m1);
    end
    
    for m2=1:num_itval
        tmp_start=t_start(m2);
        tmp_end=t_end(m2);
        ana_lfp=tmp_lfp(no_electrode,tmp_start:tmp_end);
        
        if err_flag==1
            [freq_power freq freq_err]=mtspectrumc(ana_lfp',params);
            data.freq_err=freq_err;
        else
            [freq_power freq]=mtspectrumc(ana_lfp',params);
        end
        
        data.freq_power=freq_power;
        data.freq=freq;
        spec_data{m1,m2}=data;
    end
end
%get the spectrum
spec_info.trial=use_trial;
spec_info.electrode=no_electrode;
spec_info.t_start=t_start;
spec_info.t_end=t_end;



















