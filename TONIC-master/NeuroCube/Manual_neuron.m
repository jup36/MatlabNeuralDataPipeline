function [x,y,z] = Manual_neuron(d)

r = d;
theta = pi*rand(1,1);
phi = 2*pi*rand(1,1);

x = round(r.*sin(theta).*cos(phi));
y = round(r.*sin(theta).*sin(phi));
z = round(r.*cos(theta));