%% align Xpos and Ypos to the stmLaser, stmLaserReach, stmLaserNoReach
% center yPosmm and yPosmm 
baseX = nanmean(nanmean(cell2mat(arrayfun(@(a) xPosmm(a-1500:a-500), ts.reachStart, 'un', 0)'))); 
baseY = nanmean(nanmean(cell2mat(arrayfun(@(a) yPosmm(a-1500:a-500), ts.reachStart, 'un', 0)'))); 

ts.stmLaserNoReach = ts.stmLaser(~ismember(ts.stmLaser, ts.stmLaserReach)); % get stmLaserNoReach

% xyStmLaser
xPosStmLaserC = arrayfun(@(a) xPosmm(a-1000:a+1999)-baseX, ts.stmLaser, 'un', 0); 
yPosStmLaserC = arrayfun(@(a) yPosmm(a-1000:a+1999)-baseY, ts.stmLaser, 'un', 0); 
rez.xyStmLaser = cellfun(@(a,b) [a;b], xPosStmLaserC, yPosStmLaserC, 'un', 0); 

% xyStmLaserReach
xPosStmLaserReachC = arrayfun(@(a) xPosmm(a-1000:a+1999)-baseX, ts.stmLaserReach, 'un', 0); 
yPosStmLaserReachC = arrayfun(@(a) yPosmm(a-1000:a+1999)-baseY, ts.stmLaserReach, 'un', 0); 
rez.xyStmLaserReach = cellfun(@(a,b) [a;b], xPosStmLaserReachC, yPosStmLaserReachC, 'un', 0); 

% xyStmLaserNoReach
xPosStmLaserNoReachC = arrayfun(@(a) xPosmm(a-1000:a+1999)-baseX, ts.stmLaserNoReach, 'un', 0); 
yPosStmLaserNoReachC = arrayfun(@(a) yPosmm(a-1000:a+1999)-baseY, ts.stmLaserNoReach, 'un', 0); 
rez.xyStmLaserNoReach = cellfun(@(a,b) [a;b], xPosStmLaserNoReachC, yPosStmLaserNoReachC, 'un', 0); 

%% align STR population data to the stmLaser, stmLaserReach, stmLaserNoReach
if length(ts.stmLaser)~=size(stmLaser.unitTimeTrial,3)
   error('The number of trials do NOT match!') 
end

% All STR neurons all stmLaser trials
rez.binStmLaserHz = bin1msSpkCountMat(nanmean(stmLaser.unitTimeTrial,3),50,50).*(1000/50); 
rez.binRchControlHz = bin1msSpkCountMat(nanmean(reach.unitTimeTrial,3),50,50).*(1000/50); 

rez.binStmLaserZ = binAvg1msSpkCountMat(cell2mat(stmLaser.SpkCountMatZ),50,50); 
rez.binRchControlZ = binAvg1msSpkCountMat(cell2mat(reach.SpkCountMatZ),50,50); 

% imagesc(rez.binStmLaserZ-rez.binRchControlZ)

% MSNs
rez.msnI = nanmean(rez.binStmLaserHz(:,20:30),2)<10; 