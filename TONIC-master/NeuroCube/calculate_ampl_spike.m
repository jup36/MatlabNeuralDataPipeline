function ampl_elec = calculate_ampl_spike(k,rad,neuron,electrode,gridsize)
% Calculate the amplitude of a spike averaged along the electrode surface

grid_i=-rad:gridsize:rad;
grid_j=-rad:gridsize:rad;

for i=1:length(grid_i)
    for j=1:length(grid_j)
        electrode_aux=electrode+[grid_i(i) grid_j(j) 0];
        dist=sqrt(sum((neuron-electrode_aux).^2));
        ampl(i,j)=k./(dist.^2);
    end
end
ampl_elec=mean(mean(ampl));