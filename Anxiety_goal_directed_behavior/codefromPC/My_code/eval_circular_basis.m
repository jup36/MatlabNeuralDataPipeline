function Kxxi = eval_circular_basis(x, xi, k, T)
%% create circular spline basis
%input:
%   x: position to be evaluated, [0,1]
%   xi: knots position
%   k: determine the complexity of reproducing kernel K(x,xi), default-3
%   T: period, default-1
%output:
%   basisObj: circular basis
%This script is based on the paper: 'Spline-based nonparametric regression
%for periodic functions and its application to directional tuning of
%neurons', Cari, Valerie, & Rob 2004

%code is written by Pengcheng Zhou, zhoupc1988@gmail.com

if nargin<2
    disp('The argument should have x and xi at least!');
    return;
else
    x = reshape(x,1,[]);
    xi = unique(mod(xi,2*pi));
end
if nargin<3
    k = 3;  %defaut
end
if nargin<4
    T = 1;
else
    T = abs(T);
end

kvec = (1:k)';
Kxxi = zeros(length(x), length(xi));
for m=1:length(xi)
    tmp_xi = xi(m);
    temp = bsxfun(@times, 1./(kvec.^4), cos(T*kvec*(x-tmp_xi)));
    Kxxi(:,m) = sum(temp,1);
end
