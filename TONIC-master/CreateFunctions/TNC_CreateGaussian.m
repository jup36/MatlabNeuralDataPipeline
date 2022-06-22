function [Gaussian] = TNC_CreateGaussian(Mu,Sigma,Time,dT)
% FUNCTION DETAILS: 
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI|JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________
% INPUTS: TNC_CreateGaussian(Mu,Sigma,Time,dT)
% 
% OUTPUTS: Gaussian

t = dT:dT:Time;

Gaussian = ( 1./( sqrt(2.*pi.*Sigma.*Sigma) ) ) .* exp( -(t-Mu).^2 ./ (2.*Sigma).^2 );

integral = trapz(Gaussian);
% integral = max(Gaussian);

Gaussian = Gaussian./integral;
