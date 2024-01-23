%% psthPCA on trial-averaged CTX and STR population activity 
filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};
figSavePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/pcaPSTH'; 

varName = 'rStartToPull'; 
fileName = {};      

p = parse_input_psthPCA( '/Volumes/Beefcake/Junchol_Data/JS2p0/collectData', fileName, 'rStartToPull',...
{'binSize',50,'stepSize',50,'binSizeZ',50,'useAllUnits',false,'PCsLogic',true,'PCs',5,'expVarLogic',false,'expVarCut',80,...
'FRcut',0.5,'nanTrialCut',.1,'cmap','rb','cAxis',[-3 6], 'rmvNaNtrials', true , 'rmvNaNunits', false} );  % 

% set a gaussian kernel to be convolved with the spike train (delta function)
decay    = 1;  % gaussian std
[kernel] = TNC_CreateGaussian(decay.*15, decay, decay.*30,1); % TNC_CreateGaussian(Mu,Sigma,Time,dT)

