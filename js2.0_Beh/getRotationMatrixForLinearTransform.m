%This code gets the clockwise rotation matrix to transform the (X,Y) coordinates
% such that coordinates are aligned parallel to the mouse midline.     
% usage: new coordinates are the product of the rotation matrix and the old
% coordinates.

% get data
load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','blockByblock_behavioralData_hXY_force'),'rez') % load rez

initC = {rez(:).hXYinit}; 

% get the x,y coordinates of the initial hand positions
xy = []; 
for f = 1:length(initC) 
    for b = 1:length(initC{f}) 
        nanI = cell2mat(cellfun(@(a) sum(isnan(a))>0, initC{f}{b}, 'un', 0)); 
        xy = [xy,cell2mat(initC{f}{b}(:,~nanI))]; % collect all the initial hand positions 
    end
end

x = xy(1,:);
y = xy(2,:);

xy1 = [x(-5<y&y<5); y(-5<y&y<5)]; 

% linear fit by regression
%X = [ones(length(xy1),1) xy1(1,:)'];
b = linspace(-10,10, length([-3.5:3.5]))'\[-3.5:3.5]';
a = angleTwoVectors([100;0],[100; 100*b(1)]); 
theta = deg2rad(a);

cwRot2d = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % clockwise rotation matrix
cwRot3d = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1]; % clockwise rotation matrix

% plot to check if rotation worked accurately
newxy1 = cwRotM*xy1;
hold on; 
plot(xy1(1,:),xy1(2,:),'bo') % old coordinates
plot(newxy1(1,:),newxy1(2,:),'ro') % new coordinates
hold off; 

save(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','clockwiseRotationMatrix'),'cwRot*')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thetaInDegrees = angleTwoVectors(u,v)
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    thetaInDegrees = real(acosd(CosTheta));
end

