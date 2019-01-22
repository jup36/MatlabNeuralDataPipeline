function plot2DneuralTrajAndEventMarkers( nTrjCell, neuralTrajCmap, eventMarkers, eventMarkerCmap, lineWidth, markerSize )
%figure;
hold on;
for fd = 1:length(nTrjCell)
    plot(nTrjCell{fd}(1,:),nTrjCell{fd}(2,:),'LineWidth',lineWidth,'color',neuralTrajCmap(fd,:))
    arrowh( nTrjCell{fd}(1,:), nTrjCell{fd}(2,:), neuralTrajCmap(fd,:), 175, 10:10:90)
    for evt = 1:length(eventMarkers)
        plot(nTrjCell{fd}(1,eventMarkers(evt)),nTrjCell{fd}(2,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
    end
end

%hold off;
pbaspect([1 1 1]); grid on;
%s=inputname(1); % take the input variable name as a string
%title(strcat(s));
xlabel('Dim1'); ylabel('Dim2');
end