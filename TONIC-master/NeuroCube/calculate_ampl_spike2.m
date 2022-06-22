function ampl_elec = calculate_ampl_spike2(k,rad,neuron,electrode,shift_matrix)
% Calculate the amplitude of a spike averaged along the electrode surface

electrode_grid=bsxfun(@plus,electrode',shift_matrix);
d1=bsxfun(@minus,neuron',electrode_grid);
d2=d1.^2;
d3=sum(d2,1);
dist=sqrt(d3);
ampl=k./(dist.^2);
ampl_elec=mean(ampl);

