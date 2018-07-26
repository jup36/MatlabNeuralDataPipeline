function [a_elec,b_elec] = filt_elec(diam)
% Calculate the coefficients of the electrode filter for the given diameter

rad = diam/2;
ro=72.5*1e4;
% Za=10e6;                  % Robinson 1967
Za=20e6;                    % Massobrio 2007
Rs=ro/(2*pi*diam);
% Ce=1e-9*0.2*pi*(rad)^2;   % Robinson 1967
Ce=1e-9*1.6*pi*(rad)^2;     % Massobrio 2007
Re=1.33e12/(pi*(rad)^2);
d_elec=(Rs^2*Ce)/Za+(Rs*Ce);
b_elec(1)=(Rs*Ce)/d_elec;
b_elec(2)=(Rs/Re)/d_elec;
a_elec(1)=1;
a_elec(2)=((Rs/Za)+(Rs/Re)+(Rs^2/(Re*Za)))/d_elec;