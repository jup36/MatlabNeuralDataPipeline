
function fillError(avg, ste, xval, uv_colorSTE,uv_colorMean,i)
if ~exist('i','var')
    i = 1;
    disp 'enter i value!'
end
mean_v = avg(i,:);
ste_v = ste(i,:);
xaxis_v = xval;
color_k = uv_colorSTE(i);
cd('C:\Documents and Settings\ajb110\Desktop\AlexWork\Functions\Misc Scripts')
load('colorCodes')

if (nargin < 3) xaxis_v=[1:length(mean_v)]'; end;
if (nargin < 5) color_k=1; end; 
 if size(xaxis_v,2) > 1
     xaxis_v=xaxis_v';
 end
 if size(mean_v,2) > 1
    mean_v=mean_v';
 end
 if size(ste_v,2) > 1
    ste_v=ste_v';
 end
[xaxis_v;flipdim(xaxis_v,1)];

[mean_v-ste_v/2;flipdim(mean_v+ste_v/2,1)];

fill([xaxis_v;flipdim(xaxis_v,1)],[mean_v-ste_v;flipdim(mean_v+ste_v,1)]',color_s{color_k}/255,'LineStyle','none');
hold on;
plot(xaxis_v, mean_v,'Color', color_m{uv_colorMean(i)}/255,'linewidth',2);
hold on ;

