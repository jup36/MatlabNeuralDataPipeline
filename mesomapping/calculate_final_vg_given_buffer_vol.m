function [ target_vg ] = calculate_final_vg_given_buffer_vol(buffer_vol, varargin)
% This function calculates the target viral genomes (vg) in the injectate 
% based on the provided buffer volume.
% 
% output: 
%   target_vg: calculated target viral genomes (vg) in the injectate
%
% params:
%   buffer_vol: the formulation buffer volume added for dilution (ul) (required)
%   aav_aliquot_vol: original AAV volume provided (ul)
%   inject_vol: injection volume (ul)
%   aav_titer_ml: titer of the original AAV (vg/ml)
%
% Note: Ensure to account for an average of 60-ul loss in the dead space of 
% syringe hub and needle by adding this to the final injectate. 

% Author: Junchol Park, Sep / 2024

p = parse_input_vg(varargin);
r = p.Results;

% Calculate the total volume (AAV aliquot + buffer)
total_vol = r.aav_aliquot_vol + buffer_vol;

% Calculate the target viral genomes (vg) in the injectate
target_vg = (r.aav_titer_ml / (1000 / r.aav_aliquot_vol)) * (r.inject_vol / total_vol);

    function p = parse_input_vg(varg)
        % parse input, and extract name-value pairs
        default_aav_titer_ml = 1e+13;  % default AAV titer (vg/ml)
        default_aav_aliquot_vol = 100; % default AAV volume provided (ul)
        default_inject_vol = 100;      % default injection volume (ul)

        p = inputParser; % create parser object

        addParameter(p, 'aav_titer_ml', default_aav_titer_ml);
        addParameter(p, 'aav_aliquot_vol', default_aav_aliquot_vol);
        addParameter(p, 'inject_vol', default_inject_vol);

        parse(p, varg{:});
    end
end
