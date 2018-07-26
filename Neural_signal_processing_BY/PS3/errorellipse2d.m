function errorellipse2d( data, color )
%errorellipse2d takes a 2-dimensional data (datapoint x 2 matrix), and
% generate plots of the 2-dimensional data with error (confidence) ellipses and  
% the eigenvectors indicated in the same plot

avg = nanmean(data,1);      % get the mean of each dataset
covmat = cov(data);         % get the covariance matrix (2x2 matrix)

[V,D] = eig(covmat);    % get the eigenvector and eigenvalues, for sanity check do: covmat*V == V*D

large_eigenval = max(max(D));   % identify the larger eigenvalue  
[large_eigenvec_ind_c,~] = find(D==large_eigenval);     
large_eigenvec = V(:,large_eigenvec_ind_c);     % identify the larger eigenvector 

% identify the smaller eigenvector
if large_eigenvec_ind_c == 1
    small_eigenval = max(D(:,2));
    small_eigenvec = V(:,2);
elseif large_eigenvec_ind_c == 2
    small_eigenval = max(D(:,1));
    small_eigenvec = V(:,1);
end

% get the angle from the large eigenvector 
angle = atan2(large_eigenvec(2),large_eigenvec(1));     % tan(theta) = y/x

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if angle < 0 
    angle = angle + 2*pi;
else
end

% rotation matrix
rotmat = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];

% Get the 95% confidence interval error ellipse
sqrtchisqval = sqrt(chi2inv(0.95,2));      % the square root of the chisqure statistic that corresponds to the 95 % confidence interval with 2 degrees of freedom (two are unknown)
theta_grid = linspace(0,2*pi);  % get the 100 linearly spaced points between 0 and 2*pi

X0 = avg(1);  % origin on the horizontal axis
Y0 = avg(2);  % origin on the vertical axis
a = sqrtchisqval*sqrt(large_eigenval);    % the axis length equals 2*sigma*sqrt(chisquare)
b = sqrtchisqval*sqrt(small_eigenval);    % the axis length equals 2*sigma*sqrt(chisquare)

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

% rotate 
r_ellipse = [ellipse_x_r; ellipse_y_r]' * rotmat;

%figure;
hold on;
plot(r_ellipse(:,1) + X0, r_ellipse(:,2) + Y0,'-','Color',color)
plot(data(:,1),data(:,2),'.','MarkerSize',15,'Color',color)
plot(avg(:,1),avg(:,2),'.','MarkerSize',30,'Color',color)

% Plot the eigenvectors
quiver(X0, Y0, large_eigenvec(1)*sqrt(large_eigenval), large_eigenvec(2)*sqrt(large_eigenval), '-k', 'LineWidth',2);
quiver(X0, Y0, small_eigenvec(1)*sqrt(small_eigenval), small_eigenvec(2)*sqrt(small_eigenval), '-k', 'LineWidth',2);

end

