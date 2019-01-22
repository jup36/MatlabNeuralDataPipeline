function plot2DneuralTrajAndEventMarkersCmap( neuralTrajCell, nTjCmap, eventMarkers, eventCmap, lineWidth, markerSize )
%This function takes 1-by-folds cell containing neural population
% trajectories rank-ordered by a movement kinematic variable, divided
% into a number of folds, and then trial-averaged within each fold.

% get colors
neuralTrajCmap  = TNC_CreateRBColormapJP(length(neuralTrajCell), nTjCmap); % get the neural trajectory colormap
eventMarkerCmap = TNC_CreateRBColormapJP(length(eventMarkers), eventCmap); % get the event marker colormap

if unique(cellfun(@(x) size(x,1), neuralTrajCell))==2 
elseif unique(cellfun(@(x) size(x,1), neuralTrajCell))>=3 
    neuralTrajCell = cellfun(@(x) x(1:2,:), neuralTrajCell, 'UniformOutput', false); 
end

figure;
hold on;
for fd = 1:length(neuralTrajCell)
    plot(neuralTrajCell{fd}(1,:),neuralTrajCell{fd}(2,:),'LineWidth',lineWidth,'color',neuralTrajCmap(fd,:))
    arrowh(neuralTrajCell{fd}(1,:),neuralTrajCell{fd}(2,:),neuralTrajCmap(fd,:),150,[20, 40, 60, 80]);
     for evt = 1:length(eventMarkers)
         plot(neuralTrajCell{fd}(1,eventMarkers(evt)),neuralTrajCell{fd}(2,eventMarkers(evt)),'o','MarkerSize',markerSize, 'MarkerFaceColor',eventMarkerCmap(evt,:), 'MarkerEdgeColor',eventMarkerCmap(evt,:));
     end
end
hold off;
pbaspect([1 1 1]); grid on;
xlabel('Dim1'); ylabel('Dim2'); zlabel('Dim3')

end