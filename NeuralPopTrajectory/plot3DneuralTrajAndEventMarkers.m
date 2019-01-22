function plot3DneuralTrajAndEventMarkers( neuralTrajCell, neuralTrajCmap, eventMarkers, eventMarkerCmap, lineWidth, markerSize )
%This function takes 1-by-folds cell containing neural population
%   trajectories rank-ordered by a movement kinematic variable, divided
%   into a number of folds, and then trial-averaged within each fold.

figure;
hold on;
for fd = 1:length(neuralTrajCell)
    plot3(neuralTrajCell{fd}(1,:),neuralTrajCell{fd}(2,:),neuralTrajCell{fd}(3,:),'LineWidth',lineWidth,'color',neuralTrajCmap(fd,:))
    for evt = 1:length(eventMarkers)
        plot3(neuralTrajCell{fd}(1,eventMarkers(evt)),neuralTrajCell{fd}(2,eventMarkers(evt)),neuralTrajCell{fd}(3,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
    end
end
hold off;
pbaspect([1 1 1]); grid on;
%s=inputname(1); % take the input variable name as a string
%title(s);
xlabel('Dim1'); ylabel('Dim2'); zlabel('Dim3')
end