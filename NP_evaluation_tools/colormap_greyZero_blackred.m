function map = colormap_greyZero_blackred()
    % custom colormap for inverse cdf plot from Nick Steinmetz
    map = hot(1000);
    map(1,:) = [0.4,0.4,0.4];
end