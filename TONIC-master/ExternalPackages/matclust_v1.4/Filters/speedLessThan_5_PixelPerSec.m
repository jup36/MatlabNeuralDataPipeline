    

function out = speedLessThan_5_PixelPerSec()


global clustdata;
global clustattrib;

currdir = pwd;
cd(clustattrib.currentfilepath);
cd ..
posfiles = dir('*.p');

[time, pos] = readposfiles(posfiles);
cd(currdir);
%fill in all times when th position was 0 with the nearest nonzero value
pos(:,1) = vectorfill(pos(:,1), 0,0);
pos(:,2) = vectorfill(pos(:,2), 0,0);
   
%smooth the position data
posfilt = gaussian(30*0.5, 60);
newpos = [smoothvect(pos(:,1), posfilt) smoothvect(pos(:,2), posfilt)];

%calculate velocity
vel = dist(newpos(1:end-1,:), newpos(2:end,:)) ./ (time(2:end) - time(1:end-1));
vel(end+1) = vel(end);

%lookup the velocities for all the spike times
spiketimes = clustdata.params(:,1) / 10000;
index =   lookup(spiketimes,time);
vel = vel(index);
out = (vel < 5);




function [times, pos] = readposfiles(posfiles)
times = [];
pos = [];
for i = 1:length(posfiles)
   posfilename = posfiles(i).name;
   [tmptimes, tmppos] = getposinfo(posfilename);
   times = [times; tmptimes];
   pos = [pos; tmppos];
end
times = double(times);
pos = double(pos);
[times, sortindex] = unique(times);
pos = pos(sortindex,:);
pos = pos(:,1:2);
times = times/10000;

