function ind=mea_select(fire_rate,direction, bin)
if ~exist('fire_rate','var')
    disp('no input firing rate! By default choose fire_rate in base workspace');
    try 
        fire_rate=evalin('base','fire_rate');
    catch
        disp('no firing rate');
        return;
    end
end

if ~exist('direction','var')
    disp('no input direction! By default choose direction in base workspace');
    try
        direction=evalin('base','direction');
    catch
        disp('no direction');
    end
end

if ~exist('bin','var')
    disp('By default, bin is set to 1');
    bin=1;
end

[temp1 ~]=size(fire_rate);
ind=zeros(1,temp1);
select_figure=figure;
for m=1:temp1
    mea_show_tuning(fire_rate(m,:),direction,bin);
    temp=input('Is this a good tuning curve? (y/n)?     ','s');
    if (temp=='y')|(temp=='Y')
        ind(m)=1;
        disp(m)
    end
end
return;
