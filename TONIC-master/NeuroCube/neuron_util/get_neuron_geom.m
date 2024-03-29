%  Copyright (c) California Institute of Technology, 2006 -- All Rights Reserved%  Royalty free license granted for non-profit research and educational purposes.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  get_neuron_geom %  
%  This is simply a wrapper around the loading of the geometry data.  It divides
%  the data in starting and ending points and diameters and returns them in separate%  arrays.  It converts that data from micrometers to meters for the LSA calculation.
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seg_start, seg_end, start_diams, end_diams] = get_neuron_geom(cellName, session)
    fileName = make_file_name(cellName, session, 'geom');
    geometryData = load(fileName);
    % really this is the number of lines of data
    num_rows=size(geometryData,1);
    seg_start = geometryData(1:num_rows,1:3);
    seg_end = geometryData(1:num_rows,6:8);
    start_diams = geometryData(1:num_rows,4);
    end_diams = geometryData(1:num_rows,9);
    % Points from NEURON are in micro-meters, and we return the data in meters...
    seg_start = seg_start * 1e-6;
    seg_end = seg_end * 1e-6;
    start_diams = start_diams * 1e-6;
    end_diams = end_diams * 1e-6;
