function [ ] = TNC_SS_UpdateGraphs(handles)

colormap(handles.cMap);

for i=1:numel(handles.segList)
    if handles.segList(i) < 1
        handles.segList(i) = 1;
    elseif handles.segList(i) > handles.numSegs
        handles.segList(i) = handles.numSegs;
    end
end

temp_segList = handles.segList;

% IF NO DATA IS AVAILABLE, PLOT DUMMY DATAPOINTS
if isempty(handles.featureData)  % Plot dummy Cluster graphs

    plot_dummy_datapoints(handles, temp_segList);

else % Plot real Cluster graphs

    numSegs = numel(handles.featureData.seg);
    disp(['Plotting data for shank ' num2str(handles.shankNum) ...
          ' segments ' num2str(temp_segList)]);

    for i=1:handles.numGraphs % handles.numGraphs == 3

        eval(['axes(handles.axes' num2str(i) ')']);
        cla; hold off;        
        
        % disp(['Plotting data for segment ' num2str(temp_segList(i))]);

        if  temp_segList(i) > 0        && ...
            temp_segList(i) <= numSegs && ...
            numel(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id)>1
        
            clustNums = unique(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id);

            if strcmp(handles.clusteringMode, 'mixmodel')
                noiseInds = [1:length(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id)];
                
                plot_automated_clustering_results(i, noiseInds, clustNums, temp_segList, handles);
                disp(['seg=' num2str(i) ' shank=' num2str(handles.shankNum) ...
                      ' id=' num2str(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id')]);
            else
                noiseInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==0);
                plot_results_with_none_or_manual_clustering(i, noiseInds, clustNums, temp_segList, handles);
            end
        end
    end

    % Waveform graphs
    plot_waveform_graphs(temp_segList, handles);
    
    % Cross correlation plots
    if max(clustNums)>0 && handles.showXCorr==1
        plot_cross_correlations(temp_segList, handles);
    end
end % Plot real Cluster graphs

% -------------------------------------------------------------------------------

function [] = plot_dummy_datapoints(handles, temp_segList)

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

% -------------------------------------------------------------------------------

function [] = plot_automated_clustering_results(i, noiseInds, clustNums, temp_segList, handles)

    Qthr = 75;
    Alphabet = get_alphabet();
    [normal_colors, lowQ_colors] = get_colors();

    params = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params;
    Y = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).Y;
    Q = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).Q;
    num_clusters = length(unique(Y));
    Y_unique = unique(Y);
    counts = zeros(1, length(Y_unique));
    for k=1:length(Y_unique)
        counts(k) = length(find(Y == Y_unique(k)));
    end

    for c = 1:num_clusters
        my_inds      = noiseInds((Y == Alphabet(c)) & (Q >= Qthr))';
        my_inds_lowQ = noiseInds((Y == Alphabet(c)) & (Q <  Qthr))';
        if strcmp(handles.zPlotName,'none') % 2D plot
            if length(my_inds) > 0
                if size(normal_colors, 2) >= c+1 && max(Q) > 0
                    my_color  = normal_colors{int8(c+1)};
                else
                    my_color = normal_colors{1};
                end
%               disp(['cluster=' num2str(c) ' color=' num2str(my_color) ' num_points=' num2str(length(my_inds)) ' total=' num2str(length(Y))]);
                plot(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds,handles.xPlotNum),...
                     handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds,handles.yPlotNum),... 
                     'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',my_color);
                view(2);
                hold on;
            elseif length(my_inds_lowQ) > 0
%               disp(['cluster=' num2str(c) ' color=' my_color ' num_points=' num2str(length(my_inds)) ' total=' num2str(length(Y))]);
                if size(lowQ_colors, 2) >= c+1 && max(Q) > 0
                    my_color = lowQ_colors{int8(c+1)};
                else
                    my_color = lowQ_colors{1};
                end
                plot(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds_lowQ,handles.xPlotNum),...
                     handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds_lowQ,handles.yPlotNum),...
                     'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',my_color);
                view(2);
                hold on;
            else
                hold on;
            end
            if max(clustNums)>0
                for p=2:numel(clustNums)
                    currInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));
                    plot(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.xPlotNum),...
                         handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.yPlotNum),...
                         'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.cMap(clustNums(p),:));
                    view(2);
                    hold on;
                end
            else
                % disp('No cluster ids were found');
            end
        else % 3D plot
            if length(my_inds) > 0
                if size(normal_colors, 2) >= c+1 && max(Q) > 0
                    my_color  = normal_colors{int8(c+1)};
                else
                    my_color = normal_colors{1};
                end
                plot3(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds,handles.xPlotNum),...
                      handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds,handles.yPlotNum),...
                      handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds,handles.zPlotNum),...
                      'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',my_color);
                view(3);
                hold on;
            elseif length(my_inds_lowQ) > 0
%               disp(['cluster=' num2str(c) ' color=' my_color ' num_points=' num2str(length(my_inds)) ' total=' num2str(length(Y))]);
                if size(lowQ_colors, 2) >= c+1 && max(Q) > 0
                    my_color = lowQ_colors{int8(c+1)};
                else
                    my_color = lowQ_colors{1};
                end
                plot3(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds_lowQ,handles.xPlotNum),...
                      handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds_lowQ,handles.yPlotNum),...
                      handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(my_inds_lowQ,handles.zPlotNum),...
                      'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',my_color);
                view(3);
                hold on;
            else
                hold on;
            end
            if max(clustNums)>0
                for p=2:numel(clustNums)
                    currInds = find(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id==clustNums(p));
                    plot3(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.xPlotNum),...
                          handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.yPlotNum),...
                          handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).params(currInds,handles.zPlotNum),...
                          'LineStyle','none','MarkerSize',handles.markSize,'Marker',handles.markStyle,'Color',handles.cMap(clustNums(p),:));
                    view(3);
                    hold on;
                end
            else
                % disp('No cluster ids were found');
            end
        end % 2D or 3D plot
    end % look through all colors / clusters

    title(['Segment ' num2str(temp_segList(i))]);
    set(gca,'TickDir','Out','color', [0.2 0.2 0.2]); box off; grid on;
    xlim = get(gca,'Xlim');
    ylim = get(gca,'Ylim');
    zlim = get(gca,'Zlim');

    if i==1
        ylabel(handles.yPlotName);
        xlabel(handles.xPlotName);
        if ~strcmp(handles.zPlotName,'none')
            zlabel(handles.zPlotName);
        end
    else
        if strcmp(handles.zPlotName,'none')
            axis([xlim ylim]);
        else
            axis([xlim ylim zlim]);
        end
    end

% -------------------------------------------------------------------------------

function [] = plot_results_with_none_or_manual_clustering(i, noiseInds, clustNums, temp_segList, handles)

    numSegs = numel(handles.featureData.seg);

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
    xlim = get(gca,'Xlim');
    ylim = get(gca,'Ylim');
    zlim = get(gca,'Zlim');

    if i==1
        ylabel(handles.yPlotName);
        xlabel(handles.xPlotName);
        if ~strcmp(handles.zPlotName,'none')
            zlabel(handles.zPlotName);
        end
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

% -------------------------------------------------------------------------------

function Alphabet = get_alphabet()
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
    return

% -------------------------------------------------------------------------------

function [normal_colors, lowQ_colors] = get_colors()
    normal_colors = {[1. 1. 1.], ... white
                     [1. 0. 0.], ... red
                     [1. .5 0.], ... orange
                     [1. 1. 0.], ... yellow
                     [.4 .8 0.], ... green
                     [0. .4 .8], ... blue
                     [0. 0. 1.], ... electric
                     [.4 0. .8], ... violet
                     [.8 0. .8], ... pink
                     [.8 0. .4], ... magenta
                     [.5 .5 .5]}; %  gray
    lowQ_colors =   {[1. 1. 1.], ... white
                     [1. .6 .6], ... light red
                     [1. .5 0.], ... light orange
                     [1. 1. 0.], ... light yellow
                     [.4 .8 0.], ... light green
                     [0. .4 .8], ... light blue
                     [0. 0. 1.], ... light electric
                     [.4 0. .8], ... light violet
                     [.8 0. .8], ... light pink
                     [.8 0. .4], ... light magenta
                     [.5 .5 .5]}; %  light gray
    return

% -------------------------------------------------------------------------------

function [] = plot_waveform_graphs(temp_segList, handles)
    i=1;
    if strcmp(handles.clusteringMode, 'mixmodel')
        Y = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).Y;
        clustNums = zeros(1, length(Y));
        for j=1:length(clustNums)
            Alphabet = get_alphabet();
            clustNums(j) = find(Alphabet == Y(j));
            [normal_colors, lowQ_colors] = get_colors();
        end
    else
        clustNums = unique(handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).id);
    end

    if max(clustNums)>0 && handles.updateWF==1

        disp('Displaying waveforms...');
        axes(handles.axes4);
        cla; hold off;
        numSites = size(handles.sessionStruct.seg(temp_segList(i)).shank(1).wfs(1).values,1);
        numSamps = size(handles.sessionStruct.seg(temp_segList(i)).shank(1).wfs(1).values,2);

        if numSites > 0 && numSamps > 0
            iSeg   = temp_segList(i);
            iShank = handles.shankNum;
            segShank =  handles.sessionStruct.seg(iSeg).shank(iShank);
            delta_y = 0;
            for p=1:numel(clustNums)
                curr_max = max(max(segShank.wfs(p).values(:,:)));
                if  delta_y < curr_max
                    delta_y = curr_max;
                end
            end
           
            if strcmp(handles.clusteringMode, 'mixmodel')
                Y = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).Y;
                num_clusters = length(unique(Y));
                for p=1:num_clusters
                    currInds = find(clustNums == p);
                    xOff = [1:numSamps] + (((p-2) * numSamps) + 10);
                    if strcmp(handles.clusteringMode, 'mixmodel')
                        my_color = normal_colors{int8(p+1)};
                    else
                        my_color = handles.cMap(clustNums(p),:);
                    end
                    if handles.avgWaves==1
                        temp_wf = zeros(numSites,numSamps);
                        for qq=1:numel(currInds)
                            if currInds(qq)
                                temp_wf = temp_wf + segShank.wfs(currInds(qq)).values(:,:);
                            end
                        end
                        for rr=1:numSites
%                           disp(['p=' num2str(p) ' rr=' num2str(rr) ' max_temp_wf=' num2str(max(temp_wf(rr,:)./numel(currInds))) ' min_temp_wf=' num2str(min(temp_wf(rr,:)./numel(currInds)))]);
                            yOff = rr * delta_y * 3.0;
                            plot(xOff , (temp_wf(rr,:)./numel(currInds)) + yOff , ...
                                 'Color', my_color, 'LineWidth',2);
                            hold on;
                        end
                    else
                        for qq=1:numel(currInds)
                            for rr=1:numSites
                                yOff = rr * delta_y * 3.0;
                                plot(xOff,segShank.wfs(currInds(qq)).values(rr,:)+yOff,...
                                     'Color', my_color);
                                hold on;
                            end
                        end
                    end
                end
                set(handles.axes4,'xtick',[],'ytick',[],'color', [0.2 0.2 0.2]);
                axes(handles.axes4); box off;
            end
        else
            disp(' ');
            disp('...Cannot display waveforms because *_ss.mat file is not available');
        end
    end

% -------------------------------------------------------------------------------

function [] = plot_cross_correlations(temp_segList, handles)
    disp('Computing cross-correlations...');
    axes(handles.axes5);
    cla; hold off;
    numSites = size(handles.sessionStruct.seg(1).shank(1).wfs(1).values,1);

    if numSites > 0
        for p=1:numel(clustNums)
            currInds1 = find(clustNums == p);
            times1 = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).ts(currInds1);
            for qq=1:numel(clustNums)
                currInds2 = find(clustNums == qq);
                times2 = handles.featureData.seg(temp_segList(i)).shank(handles.shankNum).ts(currInds2);
                [TNCxcorr_x,TNCxcorr_y] = TNC_CrossCorrFromTimesList(times1,times2,25);
%               disp(['p=' num2str(p) ' qq=' num2str(qq) ' size(TNCxcorr_x)=' num2str( size(TNCxcorr_x)) ' size(TNCxcorr_y)=' num2str(size(TNCxcorr_y))]);
                if qq==p
                    TNCxcorr_y(26) = 0;
                end
                dim = min(numel(TNCxcorr_x), numel(TNCxcorr_y));
                % Now plot the crosscorr
                plot(TNCxcorr_x(1:dim)+ ((p-2)*55), (TNCxcorr_y(1:dim) ./ (max(TNCxcorr_y) - min(TNCxcorr_y))) - ((qq-2).*2),'Color', handles.cMap(clustNums(p),:)); hold on; %
            end
        end
        set(handles.axes5,'xtick',[],'ytick',[],'color', [0.2 0.2 0.2]);
        axes(handles.axes5); axis tight;
    else
        disp(' ');
        disp('...Cannot display correlation plots because *_ss.mat file is not available');
    end
