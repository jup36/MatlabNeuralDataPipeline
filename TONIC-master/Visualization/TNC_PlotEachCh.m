function [] = TNC_PlotEachCh(elecDims,data1,data2,data3,reMap);
% FUNCTION DETAILS: Plotting of up to three data segments from each channel of a continuous recording in a single matrix-style plot.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

remSize = size(reMap);
dataDims = size(data1);

figure(1);
clf;

for index = 1:elecDims(1,1)*elecDims(1,2)
 subplot(elecDims(1,1),elecDims(1,2),index);
 hold on;
 
 if (max(remSize)>1)
    dataIndex = reMap(index);
 else
    dataIndex = index;
 end
 
 plot(1:dataDims(1,2),data1(dataIndex,:),'k');
 plot(1:dataDims(1,2),data2(dataIndex,:),'b');
 plot(1:dataDims(1,2),data3(dataIndex,:),'r');
 axis([0 dataDims(1,2) -20 20]);
 axis off;
 drawnow;
end