function [] = TNC_TrigPopPlotter(phys,popVec,evInds,window,fNum)

figure(fNum); clf;

wfLngth = size(phys.unit(1).wf,2);
yScale  = 1200;
numEv   = numel(evInds);
nPp     = 30;

for i=1:1:ceil(numEv./30)
    
    figure(fNum);
    clf;
        
    % display raster plot arranged by unit
 
    for p = 1:nPp
    
        cWinL = evInds((i.*nPp)+p)   - window(1);
        cWinR = evInds((i.*nPp)+p)   + window(2);

        subplot(2,nPp./2,p);
        
        for j=1:phys.numUnits

            plot([-window(1) window(2)],[yScale.*j yScale.*j],'Color',[0.5 0.5 0.5]);
            if (i.*nPp)+p-1>0
                title(evInds((i.*nPp)+p)-evInds((i.*nPp)+p-1));
            else
                title(evInds((i.*nPp)+p));                
            end
            hold on;
            validStamps = find(phys.unit(j).ts>(cWinL.*30) & phys.unit(j).ts<(cWinR.*30) );
            if numel(validStamps)>0
                for k=1:numel(validStamps)
                    plot(  ceil(phys.unit(j).ts(validStamps(k))./30) - cWinL .* [1 1],  [min(phys.unit(j).wf(validStamps(k),:)) max(phys.unit(j).wf(validStamps(k),:))] + (yScale.*j), 'Color', [j./phys.numUnits 0.67-(0.67*(j./phys.numUnits)) 1-(j./phys.numUnits)] );
                end
            end
        end
       
        plot(-window(1):window(2), popVec.proj(1,cWinL:cWinR).*yScale,'k');
        axis([-window(1) window(2) -yScale.*4 (phys.numUnits+1).*yScale]);
        axis off;

    end
    
    drawnow;
    pause();
    
end
