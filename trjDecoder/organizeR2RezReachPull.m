% filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};

%filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
%            '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'};

function [ctxRchR2, strRchR2, ctxstrRchR2, ctxPullR2, strPullR2, ctxstrPullR2] = organizeR2RezReachPull(filePath)

for f = 1:length(filePath) 
    %% load files
    rDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_reach_*')); 
    rS = load(fullfile(rDir.folder, rDir.name)); 
    pDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_pull_*')); 
    pS = load(fullfile(pDir.folder, pDir.name)); 
    % organize r2    
    for k = 1:length(rS.r2Rez.ctx)
        switch k
            case 1
                % reach position X
                ctxRchR2.posXyz{f,1}(1,1) = rS.r2Rez.ctx{k}.all; 
                strRchR2.posXyz{f,1}(1,1) = rS.r2Rez.str{k}.all; 
                ctxstrRchR2.posXyz{f,1}(1,1) = rS.r2Rez.ctxstr{k}.all; 
                % pull position X
                ctxPullR2.posXyz{f,1}(1,1) = pS.r2Rez.ctx{k}.all; 
                strPullR2.posXyz{f,1}(1,1) = pS.r2Rez.str{k}.all; 
                ctxstrPullR2.posXyz{f,1}(1,1) = pS.r2Rez.ctxstr{k}.all; 
            case 2
                % reach position Y
                ctxRchR2.posXyz{f,1}(1,2) = rS.r2Rez.ctx{k}.all; 
                strRchR2.posXyz{f,1}(1,2) = rS.r2Rez.str{k}.all; 
                ctxstrRchR2.posXyz{f,1}(1,2) = rS.r2Rez.ctxstr{k}.all; 
                % pull position Y
                ctxPullR2.posXyz{f,1}(1,2) = pS.r2Rez.ctx{k}.all; 
                strPullR2.posXyz{f,1}(1,2) = pS.r2Rez.str{k}.all; 
                ctxstrPullR2.posXyz{f,1}(1,2) = pS.r2Rez.ctxstr{k}.all; 
            case 3
                % reach position Z
                ctxRchR2.posXyz{f,1}(1,3) = rS.r2Rez.ctx{k}.all; 
                strRchR2.posXyz{f,1}(1,3) = rS.r2Rez.str{k}.all; 
                ctxstrRchR2.posXyz{f,1}(1,3) = rS.r2Rez.ctxstr{k}.all; 
                % pull position Z
                ctxPullR2.posXyz{f,1}(1,3) = pS.r2Rez.ctx{k}.all; 
                strPullR2.posXyz{f,1}(1,3) = pS.r2Rez.str{k}.all; 
                ctxstrPullR2.posXyz{f,1}(1,3) = pS.r2Rez.ctxstr{k}.all; 
            case 4
                % reach velocity X
                ctxRchR2.velXyz{f,1}(1,1) = rS.r2Rez.ctx{k}.all; 
                strRchR2.velXyz{f,1}(1,1) = rS.r2Rez.str{k}.all; 
                ctxstrRchR2.velXyz{f,1}(1,1) = rS.r2Rez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullR2.velXyz{f,1}(1,1) = pS.r2Rez.ctx{k}.all; 
                strPullR2.velXyz{f,1}(1,1) = pS.r2Rez.str{k}.all; 
                ctxstrPullR2.velXyz{f,1}(1,1) = pS.r2Rez.ctxstr{k}.all; 
            case 5
                % reach velocity X
                ctxRchR2.velXyz{f,1}(1,2) = rS.r2Rez.ctx{k}.all; 
                strRchR2.velXyz{f,1}(1,2) = rS.r2Rez.str{k}.all; 
                ctxstrRchR2.velXyz{f,1}(1,2) = rS.r2Rez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullR2.velXyz{f,1}(1,2) = pS.r2Rez.ctx{k}.all; 
                strPullR2.velXyz{f,1}(1,2) = pS.r2Rez.str{k}.all; 
                ctxstrPullR2.velXyz{f,1}(1,2) = pS.r2Rez.ctxstr{k}.all; 
            case 6
                % reach velocity X
                ctxRchR2.velXyz{f,1}(1,3) = rS.r2Rez.ctx{k}.all; 
                strRchR2.velXyz{f,1}(1,3) = rS.r2Rez.str{k}.all; 
                ctxstrRchR2.velXyz{f,1}(1,3) = rS.r2Rez.ctxstr{k}.all; 
                % pull velocity X
                ctxPullR2.velXyz{f,1}(1,3) = pS.r2Rez.ctx{k}.all; 
                strPullR2.velXyz{f,1}(1,3) = pS.r2Rez.str{k}.all; 
                ctxstrPullR2.velXyz{f,1}(1,3) = pS.r2Rez.ctxstr{k}.all; 
            case 7
                % reach position XYZ (fit together)
                ctxRchR2.posXyzT{f,1} = rS.r2Rez.ctx{k}.all; 
                ctxRchR2.posXyzToverall{f,1} = rS.r2Rez.ctx{k}.overall; 
                strRchR2.posXyzT{f,1} = rS.r2Rez.str{k}.all; 
                strRchR2.posXyzToverall{f,1} = rS.r2Rez.str{k}.overall; 
                ctxstrRchR2.posXyzT{f,1} = rS.r2Rez.ctxstr{k}.all; 
                ctxstrRchR2.posXyzToverall{f,1} = rS.r2Rez.ctxstr{k}.overall; 
                % pull position XYZ (fit together)
                ctxPullR2.posXyzT{f,1} = pS.r2Rez.ctx{k}.all; 
                ctxPullR2.posXyzToverall{f,1} = pS.r2Rez.ctx{k}.overall; 
                strPullR2.posXyzT{f,1} = pS.r2Rez.str{k}.all; 
                strPullR2.posXyzToverall{f,1} = pS.r2Rez.str{k}.overall; 
                ctxstrPullR2.posXyzT{f,1} = pS.r2Rez.ctxstr{k}.all; 
                ctxstrPullR2.posXyzToverall{f,1} = pS.r2Rez.ctxstr{k}.overall;          
            case 8
                % reach velocity XYZ (fit together)
                ctxRchR2.velXyzT{f,1} = rS.r2Rez.ctx{k}.all; 
                ctxRchR2.velXyzToverall{f,1} = rS.r2Rez.ctx{k}.overall; 
                strRchR2.velXyzT{f,1} = rS.r2Rez.str{k}.all; 
                strRchR2.velXyzToverall{f,1} = rS.r2Rez.str{k}.overall; 
                ctxstrRchR2.velXyzT{f,1} = rS.r2Rez.ctxstr{k}.all; 
                ctxstrRchR2.velXyzToverall{f,1} = rS.r2Rez.ctxstr{k}.overall; 
                % pull velocity XYZ (fit together)
                ctxPullR2.velXyzT{f,1} = pS.r2Rez.ctx{k}.all; 
                ctxPullR2.velXyzToverall{f,1} = pS.r2Rez.ctx{k}.overall; 
                strPullR2.velXyzT{f,1} = pS.r2Rez.str{k}.all; 
                strPullR2.velXyzToverall{f,1} = pS.r2Rez.str{k}.overall; 
                ctxstrPullR2.velXyzT{f,1} = pS.r2Rez.ctxstr{k}.all; 
                ctxstrPullR2.velXyzToverall{f,1} = pS.r2Rez.ctxstr{k}.overall;     
        end
    end
    % organize corr    
    fprintf('processed file# %d\n', f)
end    
  
      
end
