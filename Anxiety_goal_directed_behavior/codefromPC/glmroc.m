function [auc, tpr, fpr] = glmroc(r, rhat,flag)
% c0 = min(rhat); 
% c1 = max(rhat); 
% dc = 0.001; 
% c = c0:dc:c1; 
c = unique(rhat(r)); 
x1 = rhat(r); 
x0 = rhat(~r); 
y1 = hist(x1, c)/length(x1); 
y0 = hist(x0, c)/length(x0); 

tpr = 1-cumsum([0,y1]); 
fpr = 1-cumsum([0,y0]); 


x = diff(fpr); 
y = tpr(1:(end-1))+diff(tpr)/2; 
auc = abs(sum(x.*y)); 

if exist('flag','var')&&(flag==1)
    figure; 
    plot(fpr, tpr); 
end
    