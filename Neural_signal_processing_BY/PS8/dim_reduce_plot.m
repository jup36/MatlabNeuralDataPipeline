function dim_reduce_plot(x,z_hat,u1)
%INPUTS:
% x: high dimensional data, each row is a data point
% z_hat: low dimensional data, each row is a data point
% u1: principle direction / factor loading
[n,p] = size(x);
mu = mean(x);
% plot the data
plot(x(:,1),x(:,2),'k.')
hold on
% plot the mean
plot(mu(1),mu(2),'g.','markersize',20)
scale = 2*max(std(x))/norm(u1);     % note that the norm of an eigenvector is always 1
% plot the line defined by the first principal component
xVals = scale*[-u1(1) u1(1)]; % +/- 2 standard deviation
yVals = scale*[-u1(2) u1(2)]; % +/- 2 standard deviation
line(mu(1)+xVals,mu(2)+yVals,'color','k')
% plot the original data projected down into the subspace defined by PC1
dataPC1 = repmat(mu,n,1) + repmat(u1',n,1).*repmat(z_hat,1,p);
plot(dataPC1(:,1),dataPC1(:,2),'r.')
for i=1:n
    xVals = [x(i,1) dataPC1(i,1)];
    yVals = [x(i,2) dataPC1(i,2)];
    line(xVals,yVals,'color','r');
end
axis equal
end
