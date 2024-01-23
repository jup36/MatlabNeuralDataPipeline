
filePath_cg = {'/Volumes/Extreme SSD/js2p0/WR37_022619/Matfiles/BehVariablesJs.mat', ...  % Cg recording contra-Cg silencing, Trj checked
               '/Volumes/Extreme SSD/js2p0/WR38_052319/Matfiles/BehVariablesJs.mat', ...  % Cg recording contra-Cg silencing, Trj checked
               '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles/jsTime1k_Kinematics_VideoFiles.mat'}; % Dual recording with contra Cg silencing, Trj checked
                

filePath_delayed_cg = {'/Volumes/Extreme SSD/js2p0/WR37_022719/Matfiles/jsTime1k_Kinematics_VideoFiles.mat', ... % B Only contra-Cg delayed silencing, Trj checked
                       '/Volumes/Extreme SSD/js2p0/WR39_091019/Matfiles/jsTime1k_Kinematics_VideoFiles.mat', ... % B Only contra-Cg delayed silencing, Trj checked 
                       '/Volumes/Extreme SSD/js2p0/WR39_091119/Matfiles/jsTime1k_Kinematics_VideoFiles.mat', ... % B Only contra-Cg delayed silencing, Trj checked 
                       '/Volumes/Extreme SSD/js2p0/WR39_100319/Matfiles/jsTime1k_Kinematics_VideoFiles.mat', ... % B Only contra-Cg delayed silencing, Trj checked      
                       '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/jsTime1k_Kinematics_VideoFiles.mat'};    % Dual recording with contra Cg delayed silencing, Trj checked 

filePath_m1 = {'/Volumes/Extreme SSD/js2p0/WR38_050119/Matfiles/BehVariablesJs.mat', ...  % M1 silencing
               '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles/jsTime1k_Kinematics_VideoFiles.mat', ...  % Corticostriatal recording M1 silencing, Trj checked     
               '/Volumes/Extreme SSD/js2p0/WR40_082219/Matfiles/jsTime1k_Kinematics_VideoFiles.mat', ...
               '/Volumes/Extreme SSD/js2p0/WR45_030220/BehVariablesJs.mat'}; % M1 silencing

%% Cg stim without delay
for j = 1:length(filePath_cg)
   cgStimEffectAnalysisJsTrj(filePath_cg{j}) 
end
 
%% Cg delayed stim
for j = 1:length(filePath_delayed_cg)
   cgDelayedStimEffectAnalysisJsTrj(filePath_delayed_cg{j}) 
end

%% M1 stim
for j = 1:length(filePath_m1)
   m1StimEffectAnalysisJsTrj(filePath_m1{j}) 
end