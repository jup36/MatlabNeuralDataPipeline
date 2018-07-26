function output_coord = Rotate_coordinates(input_coord,model)
% Rotate the coordinates of each neuron to give the right orientation

xi=input_coord(1);
yi=input_coord(2);
zi=input_coord(3);
switch model
    case 1
        xo=zi;
        yo=yi;
        zo=-xi;
    case 2
        xo=-zi;
        yo=yi;
        zo=xi;
    case 3
        xo=xi;
        yo=yi;
        zo=zi;
    case 4
        xo=xi;
        yo=-zi;
        zo=yi;
    case 5
        xo=xi;
        yo=yi;
        zo=zi;
    case 6
        xo=xi;
        yo=yi;
        zo=zi;
end
output_coord(1)=xo;
output_coord(2)=yo;
output_coord(3)=zo;

return