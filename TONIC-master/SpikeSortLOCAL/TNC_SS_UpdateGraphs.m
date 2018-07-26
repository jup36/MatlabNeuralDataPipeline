function [] = TNC_SS_UpdateGraphs(handles)

colormap(handles.cMap);

for i=1:numel(handles.segList)
    if handles.segList(i) < 1
        handles.segList(i) = 1;
    elseif handles.segList(i) > handles.numSegs
        handles.segList(i) = handles.numSegs;
    end            
end
    
temp_segList = handles.segList;

if isempty(handles.featureData)

    for i=1:handles.numGraphs

        eval(['axes(handles.axes' num2str(i) ')']);
        cla;

        if strcmp(handles.zPlotName,'none')
            plot(randn(1,300)+[ones(1,100) 5.0.*ones(1,100) 2.0.*ones(1,100)], ...
                 rand(1,300)+[1.0.*ones(1,100) 2.0.*ones(1,100) 5.*ones(1,100)],...
                 'LineStyle','none','MarkerSize',handles.markSize,'Marker',...
                 handles.markStyle,'Color',handles.allFeatColor); 
        else
            plot3(randn(1,300)+[ones(1,100) 5.0.*ones(1,100) 2.0.*ones(1,100)], ...
                 rand(1,300)+[1.0.*ones(1,100) 2.0.*ones(1,100) 5.*ones(1,100)], ...
                 randn(1,300)+[5.*ones(1,100) -5.*ones(1,100) ones(1,100)],...
                 'LineStyle','none','MarkerSize',handles.markSize,'Marker',...
                 handles.markStyle,'Color',handles.allFeatColor); 
        end
        
        title(['Segment ' num2str(temp_segList(i))]);
        set(gca,'TickDir','out'); box off; grid on;
        
    end

else

    disp(['Plotting data for shank ' num2str(handles.shankNum)]);
    numSegs = numel(handles.featureData.seg);
    
    for i=1:handles.numGraphs

        eval(['axes(handles.axes' num2str(i) ')']);
        cla; hold off;        
        
        % disp(['Plotting data for segment ' num2str(temp_segList(i))]);

        if temp_segList(i) > 0 && ...
           temp_segList(i) <= numSegs && ...
           numel(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id)>0
        
            noiseInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==0);
            clustNums = unique(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id);

            if strcmp(handles.zPlotName,'none')
                
                plot(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(noiseInds,handles.xPlotNum),...
                     handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(noiseInds,handles.yPlotNum),...
                     'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.allFeatColor); 
                     hold on;

                if max(clustNums)>0
                    for p=2:numel(clustNums)                        
                        currInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));
                        plot(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.xPlotNum),...
                             handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.yPlotNum),...
                             'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.cMap(clustNums(p),:)); 
                             hold on;
                    end
                else
                    % disp('No cluster ids were found');
                end

            else
                %plot3(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(:,handles.xPlotNum),handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(:,handles.yPlotNum),handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(:,handles.zPlotNum),'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.allFeatColor); 
                plot3(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(noiseInds,handles.xPlotNum),...
                      handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(noiseInds,handles.yPlotNum),...
                      handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(noiseInds,handles.zPlotNum),...
                      'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.allFeatColor); 
                      hold on;

                if max(clustNums)>0
                    for p=2:numel(clustNums)                        
                        currInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));
                        plot3(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.xPlotNum),...
                              handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.yPlotNum),...
                              handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.zPlotNum),...
                              'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.cMap(clustNums(p),:)); 
                              hold on;
                    end
                else
                    % disp('No cluster ids were found');
                end
                
            end

            title(['Segment ' num2str(temp_segList(i))]);
            set(gca,'TickDir','Out','color', [0.2 0.2 0.2]); box off; grid on;

            if i==1
                ylabel(handles.yPlotName);
                xlabel(handles.xPlotName);
                if ~strcmp(handles.zPlotName,'none')
                    zlabel(handles.zPlotName);
                end
                xlim = get(gca,'Xlim');
                ylim = get(gca,'Ylim');
                zlim = get(gca,'Zlim');
            else
                if strcmp(handles.zPlotName,'none')
                    axis([xlim ylim]);
                else
                    axis([xlim ylim zlim]);
                end
            end

            if max(clustNums)>0
                for p=2:numel(clustNums)

                    % add labels
                    if handles.lblTrue(i)==1 && i==1
                        spcX = (xlim(2)-xlim(1))./50;
                        spcY = (ylim(2)-ylim(1))./25;
                        text(xlim(1)+spcX,ylim(2)-(spcY.*(p-1)),num2str(clustNums(p)),'BackgroundColor',handles.cMap(clustNums(p),:));
                    end

                    % plot centers
                    if handles.cntTrue(i) == 1
                        [handles]   = TNC_SS_UpdateClusterCenters(handles);
                        plot(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).cnt(clustNums(p),handles.xPlotNum),...
                             handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).cnt(clustNums(p),handles.yPlotNum),...
                             'LineStyle','none','MarkerSize',20,'Marker','+','Color',[1 1 1]);
                    end

                    % plot boundaries
                    if handles.conTrue(i) == 1
                        [handles]   = TNC_SS_UpdateClusterCenters(handles);
                        [eX,eY]     = TNC_SS_UpdateClusterBoundaries(handles,temp_segList(i),clustNums(p),handles.boundMethod);
                        plot(eX,eY,'-','Color',handles.cMap(clustNums(p),:)./2,'LineWidth',2);
                    end
                end
            end
        end
    end
    
    % Waveform graphs
    i=1;
    clustNums = unique(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id);

    if max(clustNums)>0 && handles.updateWF==1

        disp('Displaying waveforms...');
        axes(handles.axes4);
        cla; hold off;
        numSites = size(handles.sessionStruct.seg(temp_segList(i)).shank(1).wfs(1).values,1);
        numSamps = size(handles.sessionStruct.seg(temp_segList(i)).shank(1).wfs(1).values,2);
        
        if isfield(handles.sessionStruct, 'resolution')
            resolution = handles.sessionStruct.resolution;
        else
            resolution = 1;
        end
        
        for p=2:numel(clustNums)                        
            currInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));
            if handles.overlay==1
                xOff = [1:numSamps] ;
            else
                xOff = [1:numSamps] + (((p-2) * numSamps) + 10);
            end
            
            if handles.avgWaves==1
                temp_wf = zeros(numSites,numSamps,'single');
                for qq=1:numel(currInds)
%                     DEBUGGING ONLY
%                     handles.shankNum;
%                     temp_segList(i);
%                     currInds(qq);
                    temp_wf = temp_wf + (single(handles.sessionStruct.seg(temp_segList(i)).shank(handles.shankNum).wfs(currInds(qq)).values(:,:)).*resolution);
                end
                for rr=1:numSites
                    yOff = rr * 500;
                    plot(xOff , (temp_wf(rr,:)./numel(currInds)) + yOff , 'Color' , handles.cMap(clustNums(p),:) , 'LineWidth' , 2); hold on; 
                end
            else
                for qq=1:numel(currInds)
                    for rr=1:numSites
                        yOff = rr * 1000;
                        plot(xOff , (single(handles.sessionStruct.seg(temp_segList(i)).shank(handles.shankNum).wfs(currInds(qq)).values(rr,:)).*resolution) + yOff ,'Color' , handles.cMap(clustNums(p),:))+(rand(1,3)); hold on;                    
                    end
                end
            end
        end
        set(handles.axes4,'xtick',[],'ytick',[],'color', [0.2 0.2 0.2]); 
        axes(handles.axes4); box off;
    end
    
    % Cross correlation plots
    if max(clustNums)>0 && handles.showXCorr==1

        axes(handles.axes5);
        cla; hold off;
        numSites = size(handles.sessionStruct.seg(1).shank(1).wfs(1).values,1);

        for p=2:numel(clustNums)                        
            currInds1 = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));
            times1 = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).ts(currInds1);
            for qq=2:numel(clustNums)
                currInds2 = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(qq));
                times2 = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).ts(currInds2);
                [TNCxcorr_x,TNCxcorr_y] = TNC_CrossCorrFromTimesList(times1,times2,25);
                if qq==p
                    TNCxcorr_y(26) = 0;
                end
                % Now plot the crosscorr
                plot(TNCxcorr_x+ ((p-2)*60), (TNCxcorr_y ./ (max(TNCxcorr_y) - min(TNCxcorr_y))) - ((qq-2).*2),'Color', handles.cMap(clustNums(p),:)); hold on; %  
            end
        end
        set(handles.axes5,'xtick',[],'ytick',[],'color', [0.2 0.2 0.2]); 
        axes(handles.axes5); axis tight;
    end
    
end

