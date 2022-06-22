%  Copyright (c) California Institute of Technology, 2006 -- All Rights Reserved
%  Royalty free license granted for non-profit research and educational purposes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  get_phi
%
%  This function caculates the potentials at a set of 3D pts, given by
%  pt_coord = [x, y, z].  The potential is produced by a set of line segments
%  with current I, lengths ds, and distances to the pt_coord (in a the reference
%  frame of each line) given by h,R.  For full explanation of the calculation
%  and notation refer to the thesis of Gary Holt, Appendix C.
%  
%  This implementation was written by Zoran Nenadic, Caltech, 12/17/2001
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = get_phi1(h,R,ds,I, sigma, pt_coord, dist, ind1maxd, ind2maxd, Par_sim)

[N n] = size(h);                    %number of segments & points

dss = ds * ones(1,n);
L = h + dss;             %calculate L

I1 = [h < 0 & L < 0];
I3 = [h > 0 & L > 0];
I2 = ~I1 & ~I3;

i1=find(I1~=0);
i2=find(I2~=0);
i3=find(I3~=0);

phi=zeros(N,n);
phi1=ones(N,n);
% N = number of segments
% n = number of grid points for electrode of finite size

% when comparing to holt thesis, note that what we call "R" here is really r^2 in Gary's notation

if (length(i1) ~= 0)
        phi(i1)= log((sqrt(h(i1).^2+R(i1))-h(i1))./(sqrt(L(i1).^2+R(i1))-L(i1)));
end

if (length(i2) ~= 0)
        phi(i2)=log(((sqrt(h(i2).^2+R(i2))-h(i2)).* (sqrt(L(i2).^2+R(i2))+L(i2)))./R(i2));
end

if (length(i3) ~= 0)
        phi(i3)= log((sqrt(L(i3).^2+R(i3))+L(i3))./(sqrt(h(i3).^2+R(i3))+h(i3)));
end

dist2 = dist*ones(N,1);

%disp(['size(phi)=', num2str(size(phi))]);
%disp(['i1=' num2str(size(i1)) ' size(i2)=' num2str(size(i2)) ' size(i3)=' num2str(size(i3)) ' size(phi)=' num2str(size(phi))]);
if ~Par_sim.point_source
    Phi  = 1/(4*pi*sigma) * (((I./ds)   *ones(1,n)).*phi);
else
    Phi = 1/(4*pi*sigma) * (((((I(ind1maxd)+I(ind2maxd))/2.*ones(N,1))./dist2)*ones(1,n)).*phi1);
%   Phi = 1/(4*pi*sigma) * (((-abs(I)./dist2)*ones(1,n)).*phi1);
end

%disp(['size(I)=' num2str(size(I)) ' size(ds)=' num2str(size(ds)) ' size(phi)=' num2str(size(phi)) ' sizesize(Phi)=' num2str(size(Phi))]);

%disp(['dist=' num2str(dist) ' min(ds)=' num2str(min(ds)) ' max(ds)=' num2str(max(ds)) ' min(phi)=' num2str(min(min(phi))) ' max(phi)=' num2str(max(max(phi))) ]);
%disp(['dist=' num2str(dist) ' mean(ds)=' num2str(mean(ds)) ' mean(phi)=' num2str(mean(mean(phi))) ' mean(ds)/mean(phi)=' num2str(mean(ds)/mean(mean(phi)))]);
%disp(['min(Phi)=' num2str(min(min(Phi))) ' max(Phi)=' num2str(max(max(Phi))) 'min(Phi1)=' num2str(min(min(Phi1))) ' max(Phi1)=' num2str(max(max(Phi1))) 'min(I)=' num2str(min(I)) ' max(I)=' num2str(max(I))]);

%disp(['mean(sum(Phi))=' num2str(mean(sum(Phi))) ' mean(sum(Phi1))=' num2str(mean(sum(Phi1)))]);
f = sum(Phi);   % 1 x n matrix

