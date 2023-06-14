function [ buffer_vol ] = formulation_buffer_volume(varargin)
%This function calculates the formulation buffer volume to be added for
% dilution that achieves the target viral genomes at the target injection  
% volume. Note that there is on average 60-ul loss in the dead space of 
% syringe hub and needle; ensure to add 60-ul in the final injectate. 

% output: 
%   buffer_vol: the formulation buffer volume 
% params:
%   aav_aliquot_vol: the original AAV volume provided (ul)
%   target_vg: the target vg in the injectate (genomes)
%   inject_vol: injection volume (ul)

% Junchol Park Feb / 2023


p = parse_input_titer(varargin);
r = p.Results; 

dil_factor = (r.aav_titer_ml/ (1000/r.aav_aliquot_vol)) / r.target_vg; 
buffer_vol = dil_factor * r.inject_vol - r.aav_aliquot_vol; 
assert(buffer_vol >= 0, 'Check the variables as the buffer volume cannot be negative!')

    function p = parse_input_titer(varg)
        % parse input, and extract name-value pairs
        default_aav_titer_ml = 1e+13;  % the original AAV titer (vg/ml)
        default_aav_aliquot_vol = 100; % the original AAV volume provided (ul)
        default_target_vg = 2e+11;  % the target vg in the injectate (genomes)
        default_inject_vol = 100; % injection volume (ul)

        p = inputParser; % create parser object

        addParameter(p, 'aav_titer_ml', default_aav_titer_ml);
        addParameter(p, 'aav_aliquot_vol', default_aav_aliquot_vol);
        addParameter(p, 'target_vg', default_target_vg);
        addParameter(p, 'inject_vol', default_inject_vol); 

        parse(p, varg{:})
    end
end