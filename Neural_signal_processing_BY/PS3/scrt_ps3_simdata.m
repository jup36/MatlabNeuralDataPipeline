
cd('C:\Users\jup36\Dropbox\NSP\PS3')

load('ps3_simdata')

%% Plot the data points in a two-dimensional space. For classes k = 1, 2, 3,
% use a red ×, green +, and blue ? for each data point, respectively. Then, set the
% axis limits of the plot to be between 0 and 20.

c1 = [];   c2 = [];   c3 = [];   

for i = 1:size(trial,1)     % number of trial (20)
    for j = 1:size(trial,2)     % number of class (3)
     switch j
         case 1
            c1 = [c1; trial(i,j).x'];
         case 2
            c2 = [c2; trial(i,j).x'];
         case 3
            c3 = [c3; trial(i,j).x'];
    end
    end
    clearvars j
end
clearvars i

%% Find the ML model parameters using results from Problem 1.
pc1 = 1/3;   % MLE for p(C1) prior probability
pc2 = 1/3;   % MLE for p(C2)
pc3 = 1/3;   % MLE for p(C3)

muc1 = (sum(c1)./size(c1,1));  % MLE for mu(C1), this should be applicable for both gaussian and poisson
muc2 = (sum(c2)./size(c2,1));  % MLE for mu(C2)
muc3 = (sum(c3)./size(c3,1));  % MLE for mu(C3)

covc1 = cov(c1);   %  var for pc1
covc2 = cov(c2);   %  var for pc1
covc3 = cov(c3);   %  var for pc1

sharedvar = pc1*covc1 + pc2*covc2 + pc3*covc3;

%% For each class, plot the ML covariance using an ellipse of the appropriate color
hold on
errorellipse2d(c1,'r')
errorellipse2d(c2,'g')
errorellipse2d(c3,'b')
axis([-5 25 -5 25])
set(gca,'TickDir','out')

%% Draw ellipse using the contour
x1 = [-5:0.2:25];
x2 = [-5:0.2:25];
[X1, X2] = meshgrid(x1,x2);

F1 = mvnpdf([X1(:) X2(:)], muc1, covc1);
F1 = reshape(F1,length(x2),length(x1));
contour(x1,x2,F1,[.0001,.0009],'r');








