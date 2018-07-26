function [winReach] = TNC_ReachVigorWindow(reachVar,reachNum,fNum)

numReaches = reachNum;

% if numReaches>10
%     p = polyfit(1:numReaches,reachVar',20);
% else
%     p = polyfit(1:numReaches,reachVar',1);
% end
% winReach.polyFit    = polyval(p,1:numReaches);

if numReaches>21
    
    winReach.sgSmth     = sgolayfilt(reachVar,3,21);
    winReach.medF       = medfilt1(reachVar,9);

else
    
    if numReaches<10
        winReach.sgSmth     = reachVar;
    else
        if mod(round(numReaches/2),2)~=0
            winReach.sgSmth     = sgolayfilt(reachVar,3,round(numReaches/2));
        else
            winReach.sgSmth     = sgolayfilt(reachVar,3,round(numReaches/2)+1);        
        end
    end
    
    winReach.medF       = medfilt1(reachVar,floor(numReaches/3));    

end

for i=1:numReaches
    if i<5
        winReach.min(i)   = (winReach.sgSmth(i) - std(reachVar(1:4)));
        winReach.max(i)   = (winReach.sgSmth(i) + std(reachVar(1:4)));
    elseif i>numReaches-4
        winReach.min(i) = (winReach.sgSmth(i) - std(reachVar(i-4:numReaches)));
        winReach.max(i) = (winReach.sgSmth(i) + std(reachVar(i-4:numReaches)));
    else
        winReach.min(i)     = winReach.sgSmth(i) - std(reachVar(i-4:i+4));
        winReach.max(i)     = winReach.sgSmth(i) + std(reachVar(i-4:i+4)); 
    end
end

% check the quality of the output
figure(fNum); clf;
patch([1:numReaches,numReaches:-1:1],[winReach.max,winReach.min(numReaches:-1:1)],[1 0.8 0.8]); hold on;
plot(1:numReaches,reachVar,'ko');
plot(1:numReaches,winReach.sgSmth,'r-','LineWidth',2);
axis tight;
