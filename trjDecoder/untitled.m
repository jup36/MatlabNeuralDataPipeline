%% plot concatenated X Y Z
ax = 1; 
timeX = 0;
figure; hold on;
for c = 1:size(s.dat.estStateCtxMean,2) % trial types
    for r = 1:20 % just to include the first block only per trial type
        if ~isempty(s.dat.state{r,c})
            ctxTrjY = s.dat.stateCtx{r,c}(ax,:); % Ctx X(1),Y(2),or Z(3) trj (left-right, horizontal hand position)
            strTrjY = s.dat.stateStr{r,c}(ax,:); % Str X(1),Y(2),or Z(3) trj
            actTrjY = s.dat.state{r,c}(ax,:); % actual X(1),Y(2),or Z(3) trj
            
            tempX = timeX+2:timeX+length(actTrjY)+1;
            timeX = timeX+length(actTrjY)+1; % update timeX
            
            if s.dat.pos1{r,c}==min([s.dat.pos1{:}]) % left
                if s.dat.trq{r,c}==min([s.dat.trq{:}]) % low-torque
                    plot(tempX,actTrjY,'k','LineWidth',1);
                    plot(tempX,ctxTrjY,'Color',[100 149 237]./255,'LineWidth',1);
                    plot(tempX,strTrjY,'Color',[50 205 50]./255,'LineWidth',1);
                elseif s.dat.trq{r,c}==max([s.dat.trq{:}]) % high-torque
                    plot(tempX,actTrjY,'k','LineWidth',2);
                    plot(tempX,ctxTrjY,'Color',[100 149 237]./255,'LineWidth',2);
                    plot(tempX,strTrjY,'Color',[50 205 50]./255,'LineWidth',2);
                end
            elseif s.dat.pos1{r,c}==max([s.dat.pos1{:}]) % right
                if s.dat.trq{r,c}==min([s.dat.trq{:}]) % low-torque
                    plot(tempX,actTrjY,'k','LineWidth',1);
                    plot(tempX,ctxTrjY,'Color',[100 149 237]./255,'LineWidth',1);
                    plot(tempX,strTrjY,'Color',[50 205 50]./255,'LineWidth',1);
                elseif s.dat.trq{r,c}==max([s.dat.trq{:}]) % high-torque
                    plot(tempX,actTrjY,'k','LineWidth',2);
                    plot(tempX,ctxTrjY,'Color',[100 149 237]./255,'LineWidth',2);
                    plot(tempX,strTrjY,'Color',[50 205 50]./255,'LineWidth',2);
                end
            end
        end
    end
end
hold off;
clearvars r c





