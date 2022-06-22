function [] = TNC_MoverSessionPlotter(phys,popVec,behavior,winSize,fNum)

% winSize is in ms

figure(fNum); clf;
stp     = round(winSize./3);
wfLngth = size(phys.unit(1).wf,2);
yScale  = 1200;

xMi = min(popVec.proj(1,:));
xMa = max(popVec.proj(1,:));
yMi = min(popVec.proj(2,:));
yMa = max(popVec.proj(2,:));
zMi = min(popVec.proj(3,:));
zMa = max(popVec.proj(3,:));
wMi = min(popVec.proj(4,:));
wMa = max(popVec.proj(4,:));

numStps = ceil(phys.maxTime ./30 ./ stp);

figure(fNum+1);
plot(1:phys.numUnits,popVec.pca.component(:,1),'ko-'); hold on;
plot(1:phys.numUnits,popVec.pca.component(:,2),'o-','Color',[1 0 0]);
plot(1:phys.numUnits,popVec.pca.component(:,3),'o-','Color',[0 0.67 1]);


for i=1:numStps
% for i=1:1
    
    figure(fNum);
    clf;
    
    cWinL = (stp.*(i-1))+1;
    cWinR = cWinL+winSize;
    
    % display raster plot arranged by unit
    
    for j=1:phys.numUnits
        subplot(5,6,[1:5 , 7:11, 13:17]);
        plot([0 winSize],[yScale.*j yScale.*j],'Color',[0.5 0.5 0.5]);
        title([num2str(i) ' out of ' num2str(numStps) ' total steps.']);
        hold on;
        validStamps = find(phys.unit(j).ts>(cWinL.*30) & phys.unit(j).ts<(cWinR.*30) );
        if numel(validStamps)>0
            for k=1:numel(validStamps)
                hold on;
                subplot(5,6,[1:5 , 7:11, 13:17]);
                plot(  ceil(phys.unit(j).ts(validStamps(k))./30) - cWinL .* [1 1],  [min(phys.unit(j).wf(validStamps(k),:)) max(phys.unit(j).wf(validStamps(k),:))] + (yScale.*j), 'Color', [j./phys.numUnits 0.67-(0.67*(j./phys.numUnits)) 1-(j./phys.numUnits)] );
                hold on;
                subplot(5,6,[6,12,18]);
                plot(  1:wfLngth,    phys.unit(j).wf(validStamps(k),:)+(yScale.*j), 'Color', [j./phys.numUnits 0.67-(0.67*(j./phys.numUnits)) 1-(j./phys.numUnits)] );
            end
        end
%         subplot(5,6,[1:5 , 7:11, 13:17]);
%         text(-10,(yScale.*j)+(yScale./5),num2str(j));
    end
    subplot(5,6,[1:5 , 7:11, 13:17]);
    plot(1:winSize, (behavior.egLFP(cWinL+1:cWinR)),'k');
    axis([-10 winSize -yScale (phys.numUnits+1).*yScale]);
    subplot(5,6,[6,12,18]);
    axis([-10 wfLngth -yScale (phys.numUnits+1).*yScale]); axis off;
    
    % display collapsed waveforms (check sorting quality over window)
    plot([0 1],[1 0]);

    colormap('jet'); 

    subplot(5,6,[19:23]);
    %title([num2str(i./popVec.totTime.*100) '%']);
    %scatter3(popVec.proj(1,cWinL:cWinR),popVec.proj(2,cWinL:cWinR),popVec.proj(3,cWinL:cWinR),20,1:winSize+1);
    %scatter3(1:winSize,  popVec.proj(2,cWinL+1:cWinR), popVec.proj(3,cWinL+1:cWinR),  5,1:winSize);
    %axis([0 winSize yMi yMa zMi zMa]); 
    %grid on; view([-10 35]);
    plot(1:winSize,  popVec.proj(1,cWinL+1:cWinR), 'Color', [0 0 0], 'LineWidth', 2); hold on;
    plot(1:winSize,  popVec.proj(2,cWinL+1:cWinR), 'Color', [1 0 0], 'LineWidth', 2);
    plot(1:winSize,  popVec.proj(3,cWinL+1:cWinR) , 'Color', [0 0.67 1], 'LineWidth', 2);
    axis([0 winSize -0.04 0.04]);

    % behavioral data in real time
    subplot(5,6,[25:29]);
    plot(1:winSize, abs(behavior.lickData(cWinL+1:cWinR)),'k');
    axis([-10 winSize min(behavior.lickData) max(behavior.lickData)]); axis off;
    
    % population vector trace for the given window
%     figure(fNum+1);


    drawnow;
    pause();
    
end
