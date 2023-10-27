function myColormap = generateColormap(color1, color2, nColors)
    % Input:
    % color1 and color2 are the two colors at the extremes in RGB format [R, G, B].
    % nColors is the total number of colors you want, including the two extremes.
    % 
    % Output:
    % myColormap is a nColorsx3 matrix representing the colormap.

    R = linspace(color1(1), color2(1), nColors);
    G = linspace(color1(2), color2(2), nColors);
    B = linspace(color1(3), color2(3), nColors);

    myColormap = [R', G', B'];

end
