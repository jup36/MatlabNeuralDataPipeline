function Update_progress(handles,i,total)
% Update the progress bar

set(handles.neurocube_figure,'CurrentAxes',handles.Progress)
progression = floor(i*100/total);
left = 100-progression;
imagesc([zeros(5,progression) ones(5,left)])
set(gca,'XTick',[])
axis image
colormap('hot')
pause(0.01)