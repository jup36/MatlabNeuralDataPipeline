function [ geometry ] = getimec3opt3geom
%This function returns the geometry of the imec3 option3 probe. 

% Order of the probe sites in the recording file
channels = 1:384;

% Site location in micrometers (x and y)
geometry = zeros(384, 2);
viHalf = 0:(384/2-1);
geometry(1:2:end,2) = viHalf * 20;
geometry(2:2:end,2) = geometry(1:2:end,2);
geometry(1:4:end,1) = 16; %0 before
geometry(2:4:end,1) = 48; %32 before
geometry(3:4:end,1) = 0;  %16 before
geometry(4:4:end,1) = 32; %48 before

% Reference sites are being excluded
ref_sites = [37 76 113 152 189 228 265 304 341 380];
channels(ref_sites) = []; 
geometry(ref_sites,:) = []; 

% Recording contact pad size in micrometers. Height x width
pad = [12 12];

% Default prm
um_per_pix = 20;
end

