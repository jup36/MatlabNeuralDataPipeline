% filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
%             '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
%             '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
%             '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
%             '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
%             '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
%             '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};

%filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
%            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'};

function [ctxRchCorr, strRchCorr, ctxstrRchCorr, ctxPullCorr, strPullCorr, ctxstrPullCorr] = organizeCorrRezReachPull(filePath)

for f = 1:length(filePath) 
    %% load files
    rDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_reach_rec_*')); 
    rS = load(fullfile(rDir.folder, rDir.name)); 
    pDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_pull_rec_*')); 
    pS = load(fullfile(pDir.folder, pDir.name)); 
    % organize corr    
    for k = 1:length(rS.corrRez.ctx)
        switch k
            case 1
                % reach position X
                ctxRchCorr.posXyz{f,1}(1,1) = rS.corrRez.ctx{k}.all; 
                strRchCorr.posXyz{f,1}(1,1) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.posXyz{f,1}(1,1) = rS.corrRez.ctxstr{k}.all; 
                % pull position X
                ctxPullCorr.posXyz{f,1}(1,1) = pS.corrRez.ctx{k}.all; 
                strPullCorr.posXyz{f,1}(1,1) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.posXyz{f,1}(1,1) = pS.corrRez.ctxstr{k}.all; 
            case 2
                % reach position Y
                ctxRchCorr.posXyz{f,1}(1,2) = rS.corrRez.ctx{k}.all; 
                strRchCorr.posXyz{f,1}(1,2) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.posXyz{f,1}(1,2) = rS.corrRez.ctxstr{k}.all; 
                % pull position Y
                ctxPullCorr.posXyz{f,1}(1,2) = pS.corrRez.ctx{k}.all; 
                strPullCorr.posXyz{f,1}(1,2) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.posXyz{f,1}(1,2) = pS.corrRez.ctxstr{k}.all; 
            case 3
                % reach position Z
                ctxRchCorr.posXyz{f,1}(1,3) = rS.corrRez.ctx{k}.all; 
                strRchCorr.posXyz{f,1}(1,3) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.posXyz{f,1}(1,3) = rS.corrRez.ctxstr{k}.all; 
                % pull position Z
                ctxPullCorr.posXyz{f,1}(1,3) = pS.corrRez.ctx{k}.all; 
                strPullCorr.posXyz{f,1}(1,3) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.posXyz{f,1}(1,3) = pS.corrRez.ctxstr{k}.all; 
            case 4
                % reach velocity X
                ctxRchCorr.velXyz{f,1}(1,1) = rS.corrRez.ctx{k}.all; 
                strRchCorr.velXyz{f,1}(1,1) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.velXyz{f,1}(1,1) = rS.corrRez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullCorr.velXyz{f,1}(1,1) = pS.corrRez.ctx{k}.all; 
                strPullCorr.velXyz{f,1}(1,1) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.velXyz{f,1}(1,1) = pS.corrRez.ctxstr{k}.all; 
            case 5
                % reach velocity X
                ctxRchCorr.velXyz{f,1}(1,2) = rS.corrRez.ctx{k}.all; 
                strRchCorr.velXyz{f,1}(1,2) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.velXyz{f,1}(1,2) = rS.corrRez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullCorr.velXyz{f,1}(1,2) = pS.corrRez.ctx{k}.all; 
                strPullCorr.velXyz{f,1}(1,2) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.velXyz{f,1}(1,2) = pS.corrRez.ctxstr{k}.all; 
            case 6
                % reach velocity X
                ctxRchCorr.velXyz{f,1}(1,3) = rS.corrRez.ctx{k}.all; 
                strRchCorr.velXyz{f,1}(1,3) = rS.corrRez.str{k}.all; 
                ctxstrRchCorr.velXyz{f,1}(1,3) = rS.corrRez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullCorr.velXyz{f,1}(1,3) = pS.corrRez.ctx{k}.all; 
                strPullCorr.velXyz{f,1}(1,3) = pS.corrRez.str{k}.all; 
                ctxstrPullCorr.velXyz{f,1}(1,3) = pS.corrRez.ctxstr{k}.all; 
            case 7
                % reach position XYZ (fit together)
                ctxRchCorr.posXyzT{f,1} = diag(rS.corrRez.ctx{k}.all)'; 
                strRchCorr.posXyzT{f,1} = diag(rS.corrRez.str{k}.all)'; 
                ctxstrRchCorr.posXyzT{f,1} = diag(rS.corrRez.ctxstr{k}.all)'; 
                % pull position XYZ (fit together)
                ctxPullCorr.posXyzT{f,1} = diag(pS.corrRez.ctx{k}.all)'; 
                strPullCorr.posXyzT{f,1} = diag(pS.corrRez.str{k}.all)'; 
                ctxstrPullCorr.posXyzT{f,1} = diag(pS.corrRez.ctxstr{k}.all)';        
            case 8
                % reach velocity XYZ (fit together)
                ctxRchCorr.velXyzT{f,1} = diag(rS.corrRez.ctx{k}.all)'; 
                strRchCorr.velXyzT{f,1} = diag(rS.corrRez.str{k}.all)';
                ctxstrRchCorr.velXyzT{f,1} = diag(rS.corrRez.ctxstr{k}.all)'; 
                % pull velocity XYZ (fit together)
                ctxPullCorr.velXyzT{f,1} = diag(pS.corrRez.ctx{k}.all)'; 
                strPullCorr.velXyzT{f,1} = diag(pS.corrRez.str{k}.all)'; 
                ctxstrPullCorr.velXyzT{f,1} = diag(pS.corrRez.ctxstr{k}.all)';  
        end
    end
    % organize corr    
    fprintf('processed file# %d\n', f)
end    
  
      
end
