filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};
figSavePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/kfDecodeCorrKinematics'; 

%% load corrRez structure and collect data
for f = 1:length(filePath)
    %% Position data
    % correlation for position XYZ
    dirP = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPos_WR*')); % note that all files are decoding allocentric positions and saved with names 'rezKFdecodeHTrjCtxStrPos_WR' (without 'hTrjF' in their filenames)
    corrPos = load(fullfile(dirP.folder, dirP.name),'corrRez');
    corrPos = corrPos.('corrRez');
    % cortex position data
    [c.rPosCtx{f,1}, tempPosCtxI] = max([diag(corrPos.rCtx_sm)'; diag(corrPos.rCtx)'],[],1);
    [c.rPosCtxHtq{f,1}, tempPosCtxHtqI] = max([diag(corrPos.rCtxHtq_sm)'; diag(corrPos.rCtxHtq)'],[],1);
    [c.rPosCtxLtq{f,1}, tempPosCtxLtqI] = max([diag(corrPos.rCtxLtq_sm)'; diag(corrPos.rCtxLtq)'],[],1);
    [c.rPosCtxLe{f,1}, tempPosCtxLeI] = max([diag(corrPos.rCtxLe_sm)'; diag(corrPos.rCtxLe)'],[],1);
    [c.rPosCtxRi{f,1}, tempPosCtxRiI] = max([diag(corrPos.rCtxRi_sm)'; diag(corrPos.rCtxRi)'],[],1);
    [c.rPosCtxRch{f,1}, tempPosCtxRchI] = max([diag(corrPos.rCtxRch_sm)'; diag(corrPos.rCtxRch)'],[],1); 
    [c.rPosCtxPul{f,1}, tempPosCtxPulI] = max([diag(corrPos.rCtxPul_sm)'; diag(corrPos.rCtxPul)'],[],1); 
    
    tempPvalCtx = [diag(corrPos.pCtx_sm)'; diag(corrPos.pCtx)'];
    tempPvalCtxHtq = [diag(corrPos.pCtxHtq_sm)'; diag(corrPos.pCtxHtq)'];
    tempPvalCtxLtq = [diag(corrPos.pCtxLtq_sm)'; diag(corrPos.pCtxLtq)'];
    tempPvalCtxLe = [diag(corrPos.pCtxLe_sm)'; diag(corrPos.pCtxLe)'];
    tempPvalCtxRi = [diag(corrPos.pCtxRi_sm)'; diag(corrPos.pCtxRi)'];
    tempPvalCtxRch = [diag(corrPos.pCtxRch_sm)'; diag(corrPos.pCtxRch)']; 
    tempPvalCtxPul = [diag(corrPos.pCtxPul_sm)'; diag(corrPos.pCtxPul)']; 
    
    for j = 1:3 % XYZ
        c.pPosCtx{f,1}(1,j) = tempPvalCtx(tempPosCtxI(1,j),j);
        c.pPosCtxHtq{f,1}(1,j) = tempPvalCtxHtq(tempPosCtxHtqI(1,j),j);
        c.pPosCtxLtq{f,1}(1,j) = tempPvalCtxLtq(tempPosCtxLtqI(1,j),j);
        c.pPosCtxLe{f,1}(1,j) = tempPvalCtxLe(tempPosCtxLeI(1,j),j);
        c.pPosCtxRi{f,1}(1,j) = tempPvalCtxRi(tempPosCtxRiI(1,j),j);
        c.pPosCtxRch{f,1}(1,j) = tempPvalCtxRch(tempPosCtxRchI(1,j),j); 
        c.pPosCtxPul{f,1}(1,j) = tempPvalCtxPul(tempPosCtxPulI(1,j),j); 
    end
    
    % striatum position data
    [c.rPosStr{f,1}, tempPosStrI] = max([diag(corrPos.rStr_sm)'; diag(corrPos.rStr)'],[],1);
    [c.rPosStrHtq{f,1}, tempPosStrHtqI] = max([diag(corrPos.rStrHtq_sm)'; diag(corrPos.rStrHtq)'],[],1);
    [c.rPosStrLtq{f,1}, tempPosStrLtqI] = max([diag(corrPos.rStrLtq_sm)'; diag(corrPos.rStrLtq)'],[],1);
    [c.rPosStrLe{f,1}, tempPosStrLeI] = max([diag(corrPos.rStrLe_sm)'; diag(corrPos.rStrLe)'],[],1);
    [c.rPosStrRi{f,1}, tempPosStrRiI] = max([diag(corrPos.rStrRi_sm)'; diag(corrPos.rStrRi)'],[],1);
    [c.rPosStrRch{f,1}, tempPosStrRchI] = max([diag(corrPos.rStrRch_sm)'; diag(corrPos.rStrRch)'],[],1); 
    [c.rPosStrPul{f,1}, tempPosStrPulI] = max([diag(corrPos.rStrPul_sm)'; diag(corrPos.rStrPul)'],[],1); 
    
    tempPvalStr = [diag(corrPos.pStr_sm)'; diag(corrPos.pStr)'];
    tempPvalStrHtq = [diag(corrPos.pStrHtq_sm)'; diag(corrPos.pStrHtq)'];
    tempPvalStrLtq = [diag(corrPos.pStrLtq_sm)'; diag(corrPos.pStrLtq)'];
    tempPvalStrLe = [diag(corrPos.pStrLe_sm)'; diag(corrPos.pStrLe)'];
    tempPvalStrRi = [diag(corrPos.pStrRi_sm)'; diag(corrPos.pStrRi)'];
    tempPvalStrRch = [diag(corrPos.pStrRch_sm)'; diag(corrPos.pStrRch)']; 
    tempPvalStrPul = [diag(corrPos.pStrPul_sm)'; diag(corrPos.pStrPul)']; 
    
    for j = 1:3 % XYZ
        c.pPosStr{f,1}(1,j) = tempPvalStr(tempPosStrI(1,j),j);
        c.pPosStrHtq{f,1}(1,j) = tempPvalStrHtq(tempPosStrHtqI(1,j),j);
        c.pPosStrLtq{f,1}(1,j) = tempPvalStrLtq(tempPosStrLtqI(1,j),j);
        c.pPosStrLe{f,1}(1,j) = tempPvalStrLe(tempPosStrLeI(1,j),j);
        c.pPosStrRi{f,1}(1,j) = tempPvalStrRi(tempPosStrRiI(1,j),j);
        c.pPosStrRch{f,1}(1,j) = tempPvalStrRch(tempPosStrRchI(1,j),j); 
        c.pPosStrPul{f,1}(1,j) = tempPvalStrPul(tempPosStrPulI(1,j),j); 
    end
    
    %% Velocity data
    % correlation for velocity XYZ
    dirP = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrVel*'));
    corrVel = load(fullfile(dirP.folder, dirP.name),'corrRez');
    corrVel = corrVel.('corrRez');
    % cortex velocity data
    [c.rVelCtx{f,1}, tempVelCtxI] = max([diag(corrVel.rCtx_sm)'; diag(corrVel.rCtx)'],[],1);
    [c.rVelCtxHtq{f,1}, tempVelCtxHtqI] = max([diag(corrVel.rCtxHtq_sm)'; diag(corrVel.rCtxHtq)'],[],1);
    [c.rVelCtxLtq{f,1}, tempVelCtxLtqI] = max([diag(corrVel.rCtxLtq_sm)'; diag(corrVel.rCtxLtq)'],[],1);
    [c.rVelCtxLe{f,1}, tempVelCtxLeI] = max([diag(corrVel.rCtxLe_sm)'; diag(corrVel.rCtxLe)'],[],1);
    [c.rVelCtxRi{f,1}, tempVelCtxRiI] = max([diag(corrVel.rCtxRi_sm)'; diag(corrVel.rCtxRi)'],[],1);
    [c.rVelCtxRch{f,1}, tempVelCtxRchI] = max([diag(corrVel.rCtxRch_sm)'; diag(corrVel.rCtxRch)'],[],1); 
    [c.rVelCtxPul{f,1}, tempVelCtxPulI] = max([diag(corrVel.rCtxPul_sm)'; diag(corrVel.rCtxPul)'],[],1); 
    
    tempPvalCtx = [diag(corrVel.pCtx_sm)'; diag(corrVel.pCtx)'];
    tempPvalCtxHtq = [diag(corrVel.pCtxHtq_sm)'; diag(corrVel.pCtxHtq)'];
    tempPvalCtxLtq = [diag(corrVel.pCtxLtq_sm)'; diag(corrVel.pCtxLtq)'];
    tempPvalCtxLe = [diag(corrVel.pCtxLe_sm)'; diag(corrVel.pCtxLe)'];
    tempPvalCtxRi = [diag(corrVel.pCtxRi_sm)'; diag(corrVel.pCtxRi)'];
    tempPvalCtxRch = [diag(corrVel.pCtxRch_sm)'; diag(corrVel.pCtxRch)']; 
    tempPvalCtxPul = [diag(corrVel.pCtxPul_sm)'; diag(corrVel.pCtxPul)']; 
    
    for j = 1:3 % XYZ
        c.pVelCtx{f,1}(1,j) = tempPvalCtx(tempVelCtxI(1,j),j);
        c.pVelCtxHtq{f,1}(1,j) = tempPvalCtxHtq(tempVelCtxHtqI(1,j),j);
        c.pVelCtxLtq{f,1}(1,j) = tempPvalCtxLtq(tempVelCtxLtqI(1,j),j);
        c.pVelCtxLe{f,1}(1,j) = tempPvalCtxLe(tempVelCtxLeI(1,j),j);
        c.pVelCtxRi{f,1}(1,j) = tempPvalCtxRi(tempVelCtxRiI(1,j),j);
        c.pVelCtxRch{f,1}(1,j) = tempPvalCtxRch(tempVelCtxRchI(1,j),j); 
        c.pVelCtxPul{f,1}(1,j) = tempPvalCtxPul(tempVelCtxPulI(1,j),j); 
    end
    
    % striatum velocity data
    [c.rVelStr{f,1}, tempVelStrI] = max([diag(corrVel.rStr_sm)'; diag(corrVel.rStr)'],[],1);
    [c.rVelStrHtq{f,1}, tempVelStrHtqI] = max([diag(corrVel.rStrHtq_sm)'; diag(corrVel.rStrHtq)'],[],1);
    [c.rVelStrLtq{f,1}, tempVelStrLtqI] = max([diag(corrVel.rStrLtq_sm)'; diag(corrVel.rStrLtq)'],[],1);
    [c.rVelStrLe{f,1}, tempVelStrLeI] = max([diag(corrVel.rStrLe_sm)'; diag(corrVel.rStrLe)'],[],1);
    [c.rVelStrRi{f,1}, tempVelStrRiI] = max([diag(corrVel.rStrRi_sm)'; diag(corrVel.rStrRi)'],[],1);
    [c.rVelStrRch{f,1}, tempVelStrRchI] = max([diag(corrVel.rStrRch_sm)'; diag(corrVel.rStrRch)'],[],1); 
    [c.rVelStrPul{f,1}, tempVelStrPulI] = max([diag(corrVel.rStrPul_sm)'; diag(corrVel.rStrPul)'],[],1); 
    
    tempPvalStr = [diag(corrVel.pStr_sm)'; diag(corrVel.pStr)'];
    tempPvalStrHtq = [diag(corrVel.pStrHtq_sm)'; diag(corrVel.pStrHtq)'];
    tempPvalStrLtq = [diag(corrVel.pStrLtq_sm)'; diag(corrVel.pStrLtq)'];
    tempPvalStrLe = [diag(corrVel.pStrLe_sm)'; diag(corrVel.pStrLe)'];
    tempPvalStrRi = [diag(corrVel.pStrRi_sm)'; diag(corrVel.pStrRi)'];
    tempPvalStrRch = [diag(corrVel.pStrRch_sm)'; diag(corrVel.pStrRch)']; 
    tempPvalStrPul = [diag(corrVel.pStrPul_sm)'; diag(corrVel.pStrPul)']; 
    
    for j = 1:3 % XYZ
        c.pVelStr{f,1}(1,j) = tempPvalStr(tempVelStrI(1,j),j);
        c.pVelStrHtq{f,1}(1,j) = tempPvalStrHtq(tempVelStrHtqI(1,j),j);
        c.pVelStrLtq{f,1}(1,j) = tempPvalStrLtq(tempVelStrLtqI(1,j),j);
        c.pVelStrLe{f,1}(1,j) = tempPvalStrLe(tempVelStrLeI(1,j),j);
        c.pVelStrRi{f,1}(1,j) = tempPvalStrRi(tempVelStrRiI(1,j),j);
        c.pVelStrRch{f,1}(1,j) = tempPvalStrRch(tempVelStrRchI(1,j),j); 
        c.pVelStrPul{f,1}(1,j) = tempPvalStrPul(tempVelStrPulI(1,j),j); 
    end
end
clearvars f

%% plot correlation scores
% correlation X-Pos Ctx vs Str
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosCtx,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStr,'un',0))]); 
print(fullfile(figSavePath,'corrXposCtxStr'),'-dpdf','-bestfit','-painters')
[~,c.stats.rXPosCtxStrPval,~,c.stats.rXPosCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rPosCtx,'un',0)), cell2mat(cellfun(@(a) a(1),c.rPosStr,'un',0)));
% correlation Y-Pos Ctx vs Str
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosCtx,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStr,'un',0))]); 
print(fullfile(figSavePath,'corrYposCtxStr'),'-dpdf','-bestfit','-painters')
[~,c.stats.rYPosCtxStrPval,~,c.stats.rYPosCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rPosCtx,'un',0)), cell2mat(cellfun(@(a) a(2),c.rPosStr,'un',0)));
% correlation Z-Pos Ctx vs Str
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosCtx,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStr,'un',0))]); 
print(fullfile(figSavePath,'corrZposCtxStr'),'-dpdf','-bestfit','-painters')
[~,c.stats.rZPosCtxStrPval,~,c.stats.rZPosCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rPosCtx,'un',0)), cell2mat(cellfun(@(a) a(3),c.rPosStr,'un',0)));

% correlation X-Vel Ctx vs Str
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelCtx,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStr,'un',0))]); 
print(fullfile(figSavePath,'corrXvelCtxStr'),'-dpdf','-bestfit','-painters')
[~,c.stats.rXVelCtxStrPval,~,c.stats.rXVelCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rVelCtx,'un',0)), cell2mat(cellfun(@(a) a(1),c.rVelStr,'un',0)));
% correlation Y-Vel Ctx vs Str
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelCtx,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStr,'un',0))]); 
print(fullfile(figSavePath,'corrYvelCtxStr'),'-dpdf','-bestfit','-painters')
[~,c.stats.rYVelCtxStrPval,~,c.stats.rYVelCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rVelCtx,'un',0)), cell2mat(cellfun(@(a) a(2),c.rVelStr,'un',0)));
% correlation Z-Vel Ctx vs Str
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelCtx,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStr,'un',0))]); 
print(fullfile(figSavePath,'corrZvelCtxStr'),'-dpdf','-bestfit','-painters')
[~,c.stats.rZVelCtxStrPval,~,c.stats.rZVelCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rVelCtx,'un',0)), cell2mat(cellfun(@(a) a(3),c.rVelStr,'un',0)));

% correlation X-Pos Ctx vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosCtxHtq,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrXposCtxStr_HighTorque'),'-dpdf','-bestfit','-painters')
[~,c.stats.rXPosHtqCtxStrPval,~,c.stats.rXPosHtqCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rPosCtxHtq,'un',0)), cell2mat(cellfun(@(a) a(1),c.rPosStrHtq,'un',0)));
% correlation Y-Pos Ctx vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosCtxHtq,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStrHtq,'un',0))]); 
[~,c.stats.rYPosHtqCtxStrPval,~,c.stats.rYPosHtqCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rPosCtxHtq,'un',0)), cell2mat(cellfun(@(a) a(2),c.rPosStrHtq,'un',0)));
print(fullfile(figSavePath,'corrYposCtxStr_HighTorque'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Ctx vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosCtxHtq,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStrHtq,'un',0))]); 
[~,c.stats.rZPosHtqCtxStrPval,~,c.stats.rZPosHtqCtxStrStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rPosCtxHtq,'un',0)), cell2mat(cellfun(@(a) a(3),c.rPosStrHtq,'un',0)));
print(fullfile(figSavePath,'corrZposCtxStr_HighTorque'),'-dpdf','-bestfit','-painters')

% correlation X-Pos Ctx vs Str Low Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosCtxLtq,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStrLtq,'un',0))]); 
print(fullfile(figSavePath,'corrXposCtxStr_LowTorque'),'-dpdf','-bestfit','-painters')
% correlation Y-Pos Ctx vs Str Low Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosCtxLtq,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStrLtq,'un',0))]); 
print(fullfile(figSavePath,'corrYposCtxStr_LowTorque'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Ctx vs Str Low Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosCtxLtq,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStrLtq,'un',0))]); 
print(fullfile(figSavePath,'corrZposCtxStr_LowTorque'),'-dpdf','-bestfit','-painters')

% correlation X-Vel Ctx vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelCtxHtq,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrXvelCtxStr_HighTorque'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Ctx vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelCtxHtq,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrYvelCtxStr_HighTorque'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Ctx vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelCtxHtq,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrZvelCtxStr_HighTorque'),'-dpdf','-bestfit','-painters')

% correlation X-Vel Ctx vs Str Low Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelCtxLtq,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStrLtq,'un',0))]); 
print(fullfile(figSavePath,'corrXvelCtxStr_LowTorque'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Ctx vs Str Low Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelCtxLtq,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStrLtq,'un',0))]); 
print(fullfile(figSavePath,'corrYvelCtxStr_LowTorque'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Ctx vs Str Low Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelCtxLtq,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStrLtq,'un',0))]); 
print(fullfile(figSavePath,'corrZvelCtxStr_LowTorque'),'-dpdf','-bestfit','-painters')

% correlation X-Pos Str Low Torque vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosStrLtq,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrXposStr_LowVsHighTorque'),'-dpdf','-bestfit','-painters')
% correlation Y-Pos Str Low Torque vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosStrLtq,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrYposStr_LowVsHighTorque'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Str Low Torque vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosStrLtq,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrZposStr_LowVsHighTorque'),'-dpdf','-bestfit','-painters')

% correlation X-Vel Str Low Torque vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelStrLtq,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrXvelStr_LowVsHighTorque'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Str Low Torque vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelStrLtq,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrYvelStr_LowVsHighTorque'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Str Low Torque vs Str High Torque
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelStrLtq,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStrHtq,'un',0))]); 
print(fullfile(figSavePath,'corrZvelStr_LowVsHighTorque'),'-dpdf','-bestfit','-painters')

%% cortex XYZ position reach vs pull phase
% correlation X-Pos Ctx Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosCtxRch,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosCtxPul,'un',0))]); 
[~,c.stats.rXPosCtxRchPulPval,~,c.stats.rXPosCtxRchPulStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rPosCtxRch,'un',0)), cell2mat(cellfun(@(a) a(1),c.rPosCtxPul,'un',0)));
print(fullfile(figSavePath,'corrXposCtx_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Y-Pos Ctx Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosCtxRch,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosCtxPul,'un',0))]); 
[~,c.stats.rYPosCtxRchPulPval,~,c.stats.rYPosCtxRchPulStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rPosCtxRch,'un',0)), cell2mat(cellfun(@(a) a(2),c.rPosCtxPul,'un',0)));
print(fullfile(figSavePath,'corrYposCtx_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Ctx Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosCtxRch,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosCtxPul,'un',0))]); 
[~,c.stats.rZPosCtxRchPulPval,~,c.stats.rZPosCtxRchPulStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rPosCtxRch,'un',0)), cell2mat(cellfun(@(a) a(3),c.rPosCtxPul,'un',0)));
print(fullfile(figSavePath,'corrZposCtx_RchVsPul'),'-dpdf','-bestfit','-painters')

%% cortex XYZ velocity reach vs pull phase
% correlation X-Vel Ctx Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelCtxRch,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelCtxPul,'un',0))]); 
[~,c.stats.rXVelCtxRchPulPval,~,c.stats.rXVelCtxRchPulStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rVelCtxRch,'un',0)), cell2mat(cellfun(@(a) a(1),c.rVelCtxPul,'un',0)));
print(fullfile(figSavePath,'corrXvelCtx_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Ctx Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelCtxRch,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelCtxPul,'un',0))]); 
[~,c.stats.rYVelCtxRchPulPval,~,c.stats.rYVelCtxRchPulStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rVelCtxRch,'un',0)), cell2mat(cellfun(@(a) a(2),c.rVelCtxPul,'un',0)));
print(fullfile(figSavePath,'corrYvelCtx_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Ctx Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelCtxRch,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelCtxPul,'un',0))]); 
[~,c.stats.rZVelCtxRchPulPval,~,c.stats.rZVelCtxRchPulStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rVelCtxRch,'un',0)), cell2mat(cellfun(@(a) a(3),c.rVelCtxPul,'un',0)));
print(fullfile(figSavePath,'corrZvelCtx_RchVsPul'),'-dpdf','-bestfit','-painters')

%% striatum XYZ position reach vs pull phase
% correlation X-Pos Str Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosStrRch,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStrPul,'un',0))]); 
[~,c.stats.rXPosStrRchPulPval,~,c.stats.rXPosStrRchPulStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rPosStrRch,'un',0)), cell2mat(cellfun(@(a) a(1),c.rPosStrPul,'un',0)));
print(fullfile(figSavePath,'corrXposStr_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Y-Pos Str Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosStrRch,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStrPul,'un',0))]); 
[~,c.stats.rYPosStrRchPulPval,~,c.stats.rYPosStrRchPulStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rPosStrRch,'un',0)), cell2mat(cellfun(@(a) a(2),c.rPosStrPul,'un',0)));
print(fullfile(figSavePath,'corrYposStr_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Str Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosStrRch,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStrPul,'un',0))]); 
[~,c.stats.rZPosStrRchPulPval,~,c.stats.rZPosStrRchPulStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rPosStrRch,'un',0)), cell2mat(cellfun(@(a) a(3),c.rPosStrPul,'un',0)));
print(fullfile(figSavePath,'corrZposStr_RchVsPul'),'-dpdf','-bestfit','-painters')

%% striatum XYZ velocity reach vs pull phase
% correlation X-Vel Str Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelStrRch,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStrPul,'un',0))]); 
[~,c.stats.rXVelStrRchPulPval,~,c.stats.rXVelStrRchPulStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rVelStrRch,'un',0)), cell2mat(cellfun(@(a) a(1),c.rVelStrPul,'un',0)));
print(fullfile(figSavePath,'corrXvelStr_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Str Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelStrRch,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStrPul,'un',0))]); 
[~,c.stats.rYVelStrRchPulPval,~,c.stats.rYVelStrRchPulStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rVelStrRch,'un',0)), cell2mat(cellfun(@(a) a(2),c.rVelStrPul,'un',0)));
print(fullfile(figSavePath,'corrYvelStr_RchVsPul'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Str Reach vs Pull
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelStrRch,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStrPul,'un',0))]); 
[~,c.stats.rZVelStrRchPulPval,~,c.stats.rZVelStrRchPulStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rVelStrRch,'un',0)), cell2mat(cellfun(@(a) a(3),c.rVelStrPul,'un',0)));
print(fullfile(figSavePath,'corrZvelStr_RchVsPul'),'-dpdf','-bestfit','-painters')

%% Cortex Vs Striatum XYZ Velocity reach period
% correlation X-Vel Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelCtxRch,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStrRch,'un',0))]); 
[~,c.stats.rXVelCtxVsStrRchPval,~,c.stats.rXVelCtxVsStrRchStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rVelCtxRch,'un',0)), cell2mat(cellfun(@(a) a(1),c.rVelStrRch,'un',0)));
print(fullfile(figSavePath,'corrXvelCtxVsStr_Rch'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelCtxRch,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStrRch,'un',0))]); 
[~,c.stats.rYVelCtxVsStrRchPval,~,c.stats.rYVelCtxVsStrRchStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rVelCtxRch,'un',0)), cell2mat(cellfun(@(a) a(2),c.rVelStrRch,'un',0)));
print(fullfile(figSavePath,'corrYvelCtxVsStr_Rch'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelCtxRch,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStrRch,'un',0))]); 
[~,c.stats.rZVelCtxVsStrRchPval,~,c.stats.rZVelCtxVsStrRchStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rVelCtxRch,'un',0)), cell2mat(cellfun(@(a) a(3),c.rVelStrRch,'un',0)));
print(fullfile(figSavePath,'corrZvelCtxVsStr_Rch'),'-dpdf','-bestfit','-painters')

%% Cortex Vs Striatum XYZ Position reach period
% correlation X-Pos Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosCtxRch,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStrRch,'un',0))]); 
[~,c.stats.rXPosCtxVsStrRchPval,~,c.stats.rXPosCtxVsStrRchStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rPosCtxRch,'un',0)), cell2mat(cellfun(@(a) a(1),c.rPosStrRch,'un',0)));
print(fullfile(figSavePath,'corrXposCtxVsStr_Rch'),'-dpdf','-bestfit','-painters')
% correlation Y-Pos Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosCtxRch,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStrRch,'un',0))]); 
[~,c.stats.rYPosCtxVsStrRchPval,~,c.stats.rYPosCtxVsStrRchStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rPosCtxRch,'un',0)), cell2mat(cellfun(@(a) a(2),c.rPosStrRch,'un',0)));
print(fullfile(figSavePath,'corrYposCtxVsStr_Rch'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosCtxRch,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStrRch,'un',0))]); 
[~,c.stats.rZPosCtxVsStrRchPval,~,c.stats.rZPosCtxVsStrRchStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rPosCtxRch,'un',0)), cell2mat(cellfun(@(a) a(3),c.rPosStrRch,'un',0)));
print(fullfile(figSavePath,'corrZposCtxVsStr_Rch'),'-dpdf','-bestfit','-painters')

%% Cortex Vs Striatum XYZ Position pull period
% correlation X-Pos Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rPosCtxPul,'un',0)) cell2mat(cellfun(@(a) a(1),c.rPosStrPul,'un',0))]); 
[~,c.stats.rXPosCtxVsStrPulPval,~,c.stats.rXPosCtxVsStrPulStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rPosCtxPul,'un',0)), cell2mat(cellfun(@(a) a(1),c.rPosStrPul,'un',0)));
print(fullfile(figSavePath,'corrXposCtxVsStr_Pul'),'-dpdf','-bestfit','-painters')
% correlation Y-Pos Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rPosCtxPul,'un',0)) cell2mat(cellfun(@(a) a(2),c.rPosStrPul,'un',0))]); 
[~,c.stats.rYPosCtxVsStrPulPval,~,c.stats.rYPosCtxVsStrPulStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rPosCtxPul,'un',0)), cell2mat(cellfun(@(a) a(2),c.rPosStrPul,'un',0)));
print(fullfile(figSavePath,'corrYposCtxVsStr_Pul'),'-dpdf','-bestfit','-painters')
% correlation Z-Pos Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rPosCtxPul,'un',0)) cell2mat(cellfun(@(a) a(3),c.rPosStrPul,'un',0))]); 
[~,c.stats.rZPosCtxVsStrPulPval,~,c.stats.rZPosCtxVsStrPulStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rPosCtxPul,'un',0)), cell2mat(cellfun(@(a) a(3),c.rPosStrPul,'un',0)));
print(fullfile(figSavePath,'corrZposCtxVsStr_Pul'),'-dpdf','-bestfit','-painters')

%% Cortex Vs Striatum XYZ Velocity pull period
% correlation X-Vel Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(1),c.rVelCtxPul,'un',0)) cell2mat(cellfun(@(a) a(1),c.rVelStrPul,'un',0))]); 
[~,c.stats.rXVelCtxVsStrPulPval,~,c.stats.rXVelCtxVsStrPulStats] = ttest(cell2mat(cellfun(@(a) a(1),c.rVelCtxPul,'un',0)), cell2mat(cellfun(@(a) a(1),c.rVelStrPul,'un',0)));
print(fullfile(figSavePath,'corrXvelCtxVsStr_Pul'),'-dpdf','-bestfit','-painters')
% correlation Y-Vel Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(2),c.rVelCtxPul,'un',0)) cell2mat(cellfun(@(a) a(2),c.rVelStrPul,'un',0))]); 
[~,c.stats.rYVelCtxVsStrPulPval,~,c.stats.rYVelCtxVsStrPulStats] = ttest(cell2mat(cellfun(@(a) a(2),c.rVelCtxPul,'un',0)), cell2mat(cellfun(@(a) a(2),c.rVelStrPul,'un',0)));
print(fullfile(figSavePath,'corrYvelCtxVsStr_Pul'),'-dpdf','-bestfit','-painters')
% correlation Z-Vel Ctx Vs Str Reach
plotCorrTwoGroups([cell2mat(cellfun(@(a) a(3),c.rVelCtxPul,'un',0)) cell2mat(cellfun(@(a) a(3),c.rVelStrPul,'un',0))]); 
[~,c.stats.rZVelCtxVsStrPulPval,~,c.stats.rZVelCtxVsStrPulStats] = ttest(cell2mat(cellfun(@(a) a(3),c.rVelCtxPul,'un',0)), cell2mat(cellfun(@(a) a(3),c.rVelStrPul,'un',0)));
print(fullfile(figSavePath,'corrZvelCtxVsStr_Pul'),'-dpdf','-bestfit','-painters')

%% helper function
function plotCorrTwoGroups(rhoMat)
% plot the max reach Amp
figure; 
randX = -.1 + (.1+.1)*rand(100,1);
hold on
for i = 1:size(rhoMat,1)
    plot([1 2]+randX(i), rhoMat(i,:),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',15)
    plot([1 2]+randX(i), rhoMat(i,:),':')
end
hold off
xlim([0.7 2.3])
%ylim(ylimit)
set(gca,'tickDir','out')
set(gca,'Ytick',0:0.05:1)
end

