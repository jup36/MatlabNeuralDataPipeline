test_spec=elec7_data;
num_trial=120;
num_dir=200;
dir=params.dir;
dir_ind=params.dir_ind;
no_electrode=1;

temp=mea_lfp_tuning(no_electrode,test_spec,1);%mean over trials. 
test_tuning=temp(:,params.dir_ind);
freq=test_spec{1,1}.freq;
figure;
image(log(test_tuning)*15);