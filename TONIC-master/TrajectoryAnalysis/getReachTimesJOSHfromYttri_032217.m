function [ reachStart reachStop reach0 pos1 pos2] = getReachTimes( positionData )
%Extracts best reach start times from .ns2 data unless BabData set to 1 for
%processed ContData files.
% Can be made to extract from other data types. Assumes 1kHz sampling
%   Assumes  X and Y position are channels 2 and 3 (line 13)

%outputs:
%  reachStart and reachStop are vectors of reach times
%  reach0 is the amplitude readout of the whole session
%  pos1 = all reach traces, aligned to start
%  pos2 = all reach traces, aligned to stop
% 
%  sgolay injects a negative dip before a reach as an artifact of the filter.
%  we can use that for reach detection, but makes for poor position traces

TbetweenR = 2000; %minimum time between reaches, in ms
maxDur=1500;

NsDATA.Data(2,:) = positionData(1,:);
NsDATA.Data(3,:) = positionData(2,:);

reach1 = sgolayfilt(sqrt((NsDATA.Data(3,:)-median(NsDATA.Data(3,:))).^2 ...
    + (NsDATA.Data(2,:)-median(NsDATA.Data(2,:))).^2),3,201); 
reach0 = moveavg(sqrt((NsDATA.Data(3,:)-median(NsDATA.Data(3,:))).^2 ...
    + (NsDATA.Data(2,:)-median(NsDATA.Data(2,:))).^2),10);

vel=[0, diff(reach1)];


scc=sort(reach1);    
 scd=sort(vel);
scc(scc<median(scc))=median(scc);
    cutoff = scc(round(length(scc)*.90)); %.85
    cutoffd = scd(round(length(scd)*.98)); % .95
    allValidSamps1=find(vel>cutoffd);
    allValidSamps1=allValidSamps1(allValidSamps1<(length(reach1)-300));
    allValidSamps=[];
    for jj = 1:length(allValidSamps1)
        if length(allValidSamps)==0
           if sum(reach1(allValidSamps1(jj):allValidSamps1(jj)+200)>cutoff)>20
               
               allValidSamps=[allValidSamps, allValidSamps1(jj)];
               
        end;
        else
        if (allValidSamps1(jj)-allValidSamps1(jj-1))>100
        if (allValidSamps1(jj)-allValidSamps(end)) >TbetweenR
            if sum(reach1(allValidSamps1(jj):allValidSamps1(jj)+200)>cutoff)>20
                
               allValidSamps=[allValidSamps, allValidSamps1(jj)];
               
        end; end; end
        end; end
    reachStart=[];
    for i = 1:length(allValidSamps)
        ii=allValidSamps(i);
        iii=find(reach1(1:ii)<cutoff*.4,1,'last');
        if iii>10
        if ii-iii>1000  %sometimes the joystick doesn't reset to zero and induces 20 secodn long reaches.  this prevents that
            scc2=sort((reach1(iii:ii)));
            cutoff2 = scc2(round(length(scc2)*.90)); %.85 
             iii=find(reach1(1:ii)<cutoff2,1,'last');
        end
            realStartMinus10=find(vel(1:iii)<0,10,'last');
        else realStartMinus10=ii;
        end
      reachStart=[reachStart, realStartMinus10(1)];
    end
    reachStart=reachStart([1,find(diff(reachStart)>TbetweenR)+1]);
    reachStart=unique(reachStart)+10;
      % this will eliminate pseudo starts,      
      %you might not want to do eliminate them
      pos1=[];
for i =reachStart
    if i>1000
pos1=[pos1; reach0(i-200:i+1499)-reach0(i)];
    end
end
%       figure; imagesc(pos1); colorbar
%          title('start times, aligned at 200ms')
      
 
      %%
      reachStop=[];
for i = 1:length(reachStart)
        ii=reach0(reachStart(i):reachStart(i)+maxDur);
        [~,reachPeakTime]=max(ii(1:700));
        ii2=ii(reachPeakTime:end);
        sortStop=sort(ii2);
        binSize=5;  %10 may be better
        jsBinned=hist(ii2,min(ii2):binSize:max(ii2));
        [~,cutoffStop1]=max(jsBinned(1:round(length(jsBinned)/3)));
        cutoffStop=cutoffStop1*binSize+10+min(ii2);
        
           N=150; %# of ms consecutively below threshold
            t=find(ii2<cutoffStop); %i think with new consec sampling, we can be loose on this
  x = diff(t)==1;
 f = find([false,x]~=[x,false]);
 g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
  almostEnd=t(f(2*g-1))-1;
  if length(almostEnd)==0
     almostEnd=find(ii2<cutoffStop,1);
     reachStart(i);
  end
        %thisStop= find(moveavg(diff(ii2(almostEnd:end)),10)>=0,1);
        reachStop=[reachStop,  reachPeakTime+almostEnd+reachStart(i)]; %thisStop+reachPeakTime+almostEnd+reachStart(i)];
end

  pos2=[];
for i =reachStop
    if i>1500 & i < (length(reach0)-200)
        pos2=[pos2; reach0(i-1499:i+200)-reach0(i)];
    end
end
%       figure; imagesc(pos2); colorbar
%       title('stop times, aligned at 1500ms')
      
%       figure; plot(reach0);
%       hold on; plot(reachStart, reach0(reachStart),'g*');
%       hold on; plot(reachStop, reach0(reachStop),'r*');
%       title('JS position with start and stop times');
      
%%
      %  if max(diff(endThresh))==1
           % almostEnd=endThresh(end)+reachStart(i);
            
       % else
        %iicut=find(diff(find(reach0>cutoffStop))>
        %end
       % reachStop=[reachStop,find(diff(reach0(almostEnd:...
        %    almostEnd+maxDur))>=0,1)+almostEnd];
end
   % below displacement threshold and 