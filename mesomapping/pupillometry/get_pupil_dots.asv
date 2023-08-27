% insert pupil points
function dots = get_pupil_dots(tableDlc, likelihoodcut)

dots = cell(8, 1);
dots(:) = deal({NaN});


if tableDlc{1, 'Lpupil_likelihood'} >= likelihoodcut
    dots{1, 1} = [tableDlc{1, 'Lpupil_x'}, tableDlc{1, 'Lpupil_y'}]; % L
end

if tableDlc{1, 'LDpupil_likelihood'} >= likelihoodcut
    dots{2, 1} = [tableDlc{1, 'LDpupil_x'}, tableDlc{1, 'LDpupil_y'}]; % LD
end

if tableDlc{1, 'Dpupil_likelihood'} >= likelihoodcut
    dots{3, 1} = [tableDlc{1, 'Dpupil_x'}, tableDlc{1, 'Dpupil_y'}]; % D
end

if tableDlc{1, 'DRpupil_likelihood'} >= likelihoodcut
    dots{4, 1} = [tableDlc{1, 'DRpupil_x'}, tableDlc{1, 'DRpupil_y'}]; % DR
end

if tableDlc{1, 'Rpupil_likelihood'} >= likelihoodcut
    dots{5, 1} = [tableDlc{1, 'Rpupil_x'}, tableDlc{1, 'Rpupil_y'}]; % R
end

if tableDlc{1, 'RVupil_likelihood'} >= likelihoodcut % RV
    dots{6, 1} = [tableDlc{1, 'RVupil_x'}, tableDlc{1, 'RVupil_y'}]; % R
end

if tableDlc{1, 'Vpupil_likelihood'} >= likelihoodcut % V
    dots{7, 1} = [tableDlc{1, 'Vpupil_x'}, tableDlc{1, 'Vpupil_y'}]; % R
end

if tableDlc{1, 'VLpupil_likelihood'} >= likelihoodcut % VL
    dots{8, 1} = [tableDlc{1, 'VLpupil_x'}, tableDlc{1, 'VLpupil_y'}]; % R
end

end
