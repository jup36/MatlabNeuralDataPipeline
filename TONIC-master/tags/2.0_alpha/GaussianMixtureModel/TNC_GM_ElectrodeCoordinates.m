function [XE, YE, ZE] = TNC_GM_ElectrodeCoordinates(probe, iE)
    if strcmp(probe, 'Buzsaki64')
        [XE, YE, ZE] = get_Buzsaki64_electrode_coordinates(iE);
    else % tetrode
        [XE, YE, ZE] = get_tetrode_electrode_coordinates(iE);
    end

%-------------------------------------------------------------------------------

function [XE, YE, ZE] = get_Buzsaki64_electrode_coordinates(iE)
    X  = [0  -7.25  7.25 -10.75  10.75 -14.5  14.5 -18];
    Y  = [22 44    66     88    110    132   154   176];
    Z  = [0   0     0      0      0      0     0     0];
    XE = X(iE);
    YE = Y(iE);
    ZE = Z(iE); 

%-------------------------------------------------------------------------------

function [XE, YE, ZE] = get_tetrode_electrode_coordinates(iE)
    X  = [ 25  25 -25 -25 ];
    Y  = [ 25 -25  25 -25 ];
    Z  = [  0   0   0   0 ];
    XE = X(iE);
    YE = Y(iE);
    ZE = Z(iE);

