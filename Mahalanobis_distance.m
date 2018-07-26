%% Mahalanobis distance
% generate correlated bivariate data
X = mvnrnd( [0;0], [1 .9; .9 1], 100 ); % random vectors from the multivariate normal distribution

% input observations
Y = [ 1 1; 1 -1; -1 1; -1 -1 ]; 

% compute the mahalanobis distance of observations in Y from the reference sample in X
d1 = mahal(Y,X); % mahalanobis distance which takes the covariance of data into account

% compute their squared euclidean distances from the mean of X for a comparison
d2 = sum((Y-repmat(mean(X),4,1)).^2, 2); % euclidean distance agnostic to the covariance of the data

% plot the observations with Y values colored according to the Mahalanobis distance 
scatter(X(:,1),X(:,2)) % scatter the reference observations
hold on
scatter(Y(:,1),Y(:,2),100,d1,'*','LineWidth',2) % scatter the input observations with the size and color of the asterisks specified
hb = colorbar;
ylabel(hb,'Mahalanobis Distance')
legend('X','Y','Location','NW')





