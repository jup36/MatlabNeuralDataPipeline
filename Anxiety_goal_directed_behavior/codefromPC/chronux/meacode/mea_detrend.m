function detrend_data=mea_detrend(ns2,params)

try 
    Fs=params.Fs;
catch
    display('You didnot not set Fs, by default, it is set to 1000Hz');
    Fs=1000;
end

try 
    mwin=params.mwin;
catch
    display('You didnot set Fs, by defaut, it is set to [0.1 0.05]');
    mwin=[0.1,0.05];
end

num_trial=length(ns2);
detrend_data=cell(1,num_trial);
for m=1:num_trial
    tmp_trial=ns2{m};
    detrend_data{m}=locdetrend(double(tmp_trial),Fs,mwin);
end
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
