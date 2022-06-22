function [] = TNC_AnimatePopVector(popVec,window,target,dims)

% get the window step size
stp = window(1);

% get the window width
width = window(2);

figure(target); clf;

xMi = min(popVec.proj(1,:));
xMa = max(popVec.proj(1,:));
yMi = min(popVec.proj(2,:));
yMa = max(popVec.proj(2,:));
zMi = min(popVec.proj(3,:));
zMa = max(popVec.proj(3,:));
wMi = min(popVec.proj(4,:));
wMa = max(popVec.proj(4,:));



for i=width+1:stp:popVec.totTime

    if dims==3
        colormap('bone'); title([num2str(i./popVec.totTime.*100) '%']);
        scatter3(popVec.proj(1,i-width:i+(width/10)),popVec.proj(2,i-width:i+(width/10)),popVec.proj(3,i-width:i+(width/10)),60,(width./10):-1:-width);
        grid on; axis([xMi xMa yMi yMa zMi zMa]);
    else
        colormap('bone');
        subplot(121);
        scatter(popVec.proj(1,i-width:i+(width/10)),popVec.proj(2,i-width:i+(width/10)),60,(width./10):-1:-width);
        grid on; axis([xMi xMa yMi yMa]); title([num2str(i./popVec.totTime.*100) '%']);
        subplot(122);
        scatter(popVec.proj(3,i-width:i+(width/10)),popVec.proj(4,i-width:i+(width/10)),60,(width./10):-1:-width);
        grid on; axis([zMi zMa wMi wMa]); title([num2str(i./popVec.totTime.*100) '%']);
    end
    drawnow;

end