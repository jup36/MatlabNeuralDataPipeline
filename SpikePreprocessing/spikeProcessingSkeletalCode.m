addpath(genpath('/Volumes/RAID2/parkj/MATLAB'))
filePath = '/Volumes/RAID2/parkj/NeuralData/ITphys/IT01_Ldms_M1_121317/Matfiles'; 

cd(filePath)

%% get behavioral time stamps
behaviorTimestamps(filePath,'numbNeuralProbe',0,...
'numbChEachProbe',64,'XposCh',33,'YposCh',37,'soleCh',3,'lickCh',5,'laserCh',7,'numbTagLasers',60,'artifactRmv',true,'reachBeforeLastReward',true)

%% spike sorting
runJrc('/Volumes/RAID2/parkj/NeuralData/ITphys/IT01_Ldms_M1_121317/Data',{'ITstim_DMSM1_record_121317_g0_t0.imec.ap.bin'})
%jrc('manual', prmFile{iFile})

%% get PSTHs (generates binSpkCount*.mat)
PSTH_rasters( filePath,'IT01_121317',4316,'probeAngle',10,'strCtx',true,'strCtxBorder',2000,'numbSiteProbe',[384],...
'psthPlotFlag',false,'reachWin',[2e3 2e3],'rewardWin',[3e3 1e3],'tagLaserWin',[5e3 5e3])

%% run PCA on PSTHs (generates pcaPSTHbinSpkCount*.mat)
pcaPSTH( filePath, 'binSpkCountCTXIT01_121317', 'reward',...
'binSize',200,'stepSize',200,'dcFactor',50,'useAllUnits',false,'PCsLogic',true,'PCs',3,'expVarLogic',false,'expVarCut',80,...
'FRcut',1,'nanTrialCut',.1,'cmap','hot','cAxis',[-3 6]) % CTX

pcaPSTH( filePath, 'binSpkCountSTRIT01_121317', 'reward',...
'binSize',200,'stepSize',200,'dcFactor',50,'useAllUnits',false,'PCsLogic',false,'PCs',3,'expVarLogic',true,'expVarCut',80,...
'FRcut',1,'nanTrialCut',.1,'cmap','hot','cAxis',[-3 6]) % STR

%% examine laser stim effects on neural activity 
[stimEctx] = stimEffect( filePath, 'binSpkCountCTXIT01_121317', 4000,...
'reachWinEdges',-2e3+1:2e3,'reachOnTime',0,'reachDur',750,'tagWinEdges',-5e3+1:5e3,'tagOnTime',0,'tagDur',500,'lowFRcut',0.5,'colorAxis',[-3 3]);
% Use IndividualUnitPlotGramm.m to plot individual unit PSTHs based on the outcome of 'stimEffect'.  

[stimEstr] = stimEffect( filePath, 'binSpkCountSTRIT01_121317', 4000,...
'reachWinEdges',-2e3+1:2e3,'reachOnTime',0,'reachDur',750,'tagWinEdges',-5e3+1:5e3,'tagOnTime',0,'tagDur',500,'lowFRcut',0.5,'colorAxis',[-3 3]);
% Use IndividualUnitPlotGramm.m to plot individual unit PSTHs based on the outcome of 'stimEffect'.  


%% examine laser stim effects on reach behavior
[stimR] = stimEffectReach(filePath); 

%% run DCA
runDCA( filePath, 'binSpkCountCTXIT01_121317pcaPSTHreward200ms', 'binSpkCountSTRIT01_121317pcaPSTHreward200ms', 100, 'IT01_121317_DCA_reward200ms' )

%% run GPFA
runGPFA( filePath )
