function [vf, vf_weighted] = vector_fields(dPrm, normalization_factor)
% dPrm: dPrm structure
load(fullfile('/Volumes/dudmanlab/junchol/js2p0/collectData','a2dColorMap.mat'),'colormap2D') % 2d colorMap for the scatter plot
numColors = size(colormap2D,1);
rangeColor = linspace(-2,2,numColors);

x = -2:.2:2;
y = 2:-.2:-2;
timeI = -1:.02:1; 

[X,Y] = meshgrid(x, y); 

dist = @(x, y) sqrt(x^2 + y^2); 

distMat = cell2mat(reshape(arrayfun(@(a, b) dist(a, b), reshape(X, [], 1), reshape(Y, [], 1), 'un', 0), length(x), length(y))); 

vf = zeros(length(x), length(y), length(timeI)); % the vector field cell

for j = 1:length(dPrm)
    trj = dPrm(j).s.trj2d;
    sigI = dPrm(j).s.trjDistZero_sigI; 
    
    if ~isnan(sum(trj(:))) && ~isempty(trj)
        for tt = 1:length(timeI) % time steps
            if sigI(tt)
                column = cell2mat(arrayfun(@map_point_X, trj(tt,1), 'un', 0));
                row = cell2mat(arrayfun(@map_point_Y, trj(tt,2), 'un', 0));
                vf(row, column, tt) = vf(row, column, tt) + 1;
            end
        end
    end
end

vf_weighted = vf.*distMat; 

vf_mean = nanmean(vf_weighted,3)./10.*normalization_factor; 

figure; hold on; 
pbaspect([1 1 1])
for rr = 1:size(vf_mean,1)
    for cc = 1:size(vf_mean,2)
        [~,cXI] = min(abs(X(rr,cc)-rangeColor));
        [~,cYI] = min(abs(Y(rr,cc)-rangeColor));
        thisColor = [colormap2D(cXI, cYI, 1), colormap2D(cXI, cYI, 2), colormap2D(cXI, cYI, 3)];
        quiver_plot([X(rr,cc), Y(rr,cc)], vf_mean(rr,cc)+0.01, thisColor);   
    end
end
xlim([-2 2])
ylim([-2 2])
set(gca, 'TickDir', 'out')
xticks(-2:1:2)
yticks(-2:1:2)
plot([-2 2], [-2 2], ':k', 'LineWidth', .5);
plot([-2 2], [2 -2], ':k', 'LineWidth', .5);
plot([-2 2], [0 0], ':k', 'LineWidth', .5);
plot([0 0], [2 -2], ':k', 'LineWidth', .5);

hold off; 





function quiver_plot(point, vec_len, color)
% https://stackoverflow.com/questions/18776172/in-matlab-how-do-i-change-the-arrow-head-style-in-quiver-plot

ref = [point(1), 0];
rad = angleTwoVecRad(point, ref);

sign_point = sign(point);

if point(1) == 0
    xd = 0;
    yd = vec_len*sign_point(2);
else
    xd = vec_len*cos(rad)*sign_point(1);
    yd = vec_len*sin(rad)*sign_point(2);
end

point2 = [point(1)-xd, point(2)-yd];
diff_point = point-point2; 

%quiver(point2(1),point2(2),diff_point(1),diff_point(2), 1, 'LineWidth', vec_len*3, 'Color', color);
L = median(sqrt(diff_point(1).^2+diff_point(2).^2))/4;
A = 22.5;
W = vec_len*3;

quiver_tri(point2(1), point2(2), diff_point(1), diff_point(2), L, A, W, color)


function quiver_tri(x,y,u,v,varargin)
%QUIVER_TRI Quiver plot with filled triangles
%   QUIVER_TRI(x,y,u,v) plots velocity vectors with filled triangle arrows.
%   Contrary to QUIVER, QUIVER_TRI does not scale vectors before plotting
%   them. Vectors will be plotted with the same spacing and magnitude as
%   they are provided. Default head arrow size is 20% of the median
%   magnitude and default head angle is 22.5 degrees.
%
%   QUIVER_TRI(x,y,u,v,headsize) plots velocity vectors specifying the
%   headsize (same units as u and v).
%
%   QUIVER_TRI(x,y,u,v,headsize,headangle) plots velocity vectors
%   specifying the head size (same units as u and v) and the head angle
%   amplitude (in degrees).
%
%   QUIVER_TRI(x,y,u,v,headsize,headangle,width) same as before, plus
%   setting of the quiver body width
%
%   QUIVER_TRI(x,y,u,v,headsize,headangle,width,col) same as before, plus
%   setting of the color
%
%   See also QUIVER.
%
%   Author: Alessandro Masullo 2015
if nargin == 0
    load wind
    x = x(:,:,5);
    y = y(:,:,5);
    u = u(:,:,5)/7;
    v = v(:,:,5)/7;
    triL = median(sqrt(u(:).^2+v(:).^2))/5;
    triA = 22.5;
    width = 1;
    col = 'k';
end
if nargin == 4
    triL = median(sqrt(u(:).^2+v(:).^2))/5;
    triA = 22.5;
    width = 1;
    col = 'k';
end
if nargin == 6
    triL = varargin{1};
    triA = varargin{2};
    width = 1;
    col = 'k';
end
if nargin == 7
    triL = varargin{1};
    triA = varargin{2};
    width = varargin{3};
    col = 'k';
end
if nargin == 8
    triL = varargin{1};
    triA = varargin{2};
    width = varargin{3};
    col = varargin{4};
end
triA = triA/180*pi;
% Coordinates of the triangle in vertex reference system
V = [0; 0];
U = [-triL*cos(triA); triL*sin(triA)];
B = [-triL*cos(triA); -triL*sin(triA)];
C = [-triL*cos(triA); 0];
t = atan2(v(:),u(:));
figure(gcf),ish = ishold; hold on
Vt = zeros(numel(u),2);
Ut = zeros(numel(u),2);
Bt = zeros(numel(u),2);
Ct = zeros(numel(u),2);
% Change coordinate system of the triangles
for i = 1:numel(u)
    % Rotate and translate the triangle to the right position
    M = [cos(t(i)) -sin(t(i)); sin(t(i)) cos(t(i))];
    T = [x(i)+u(i); y(i)+v(i)];
    
    Vt(i,:) = M*V + T;
    Ut(i,:) = M*U + T;
    Bt(i,:) = M*B + T;
    Ct(i,:) = M*C + T;
end
% Plot the arrow lines
plot([x(:)'; Ct(:,1)'],[y(:)'; Ct(:,2)'],'-','linewidth',width,'color',col)
% Plot the triangle heads
patch([Vt(:,1)'; Ut(:,1)'; Bt(:,1)'; Vt(:,1)'],[Vt(:,2)'; Ut(:,2)'; Bt(:,2)'; Vt(:,2)'],col,'EdgeColor',col)
if ~ish
    hold off
end

end


end




function thetaInRadian = angleTwoVecRad(u,v)
    cosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    thetaInDegrees = real(acosd(cosTheta));
    thetaInRadian = deg2rad(thetaInDegrees); 
end

function minI = map_point_X(point)
ref = -2:.2:2; 
[~, minI] = min(abs(point - ref));
end

function minI = map_point_Y(point)
ref = 2:-.2:-2; 
[~, minI] = min(abs(point - ref));
end



end

