function loc = locateOnArray(array, target)
    % Convert the array to a column vector for broadcasting
    array = array(:);
    target = target(:); 

    % Subtract each target value from the entire array and find the minimum absolute value
    [~, loc] = min(abs(array - target.'), [], 1);
end
