function test_read_text(input_file)
    % Example: 
    % >> test_read_text('NeuroCubeInput.txt')

    i = 1;
    fid = fopen(input_file,'r'); %# open csv file for reading
    while ~feof(fid)
        line = fgets(fid); %# read line by line
        coords(i,1:3) = sscanf(line,'%f %f %f'); %# sscanf can read only numeric data :(
        i = i+1;
    end
    fclose(fid);
    Neurons.nneurons = size(coords, 1);
    disp(['Neurons.nneurons=' num2str(Neurons.nneurons)]);
    Neurons.x_neurons = coords(:,1)';
    Neurons.y_neurons = coords(:,2)';
    Neurons.z_neurons = coords(:,3)';
    Neurons.x_neurons
    Neurons.y_neurons
    Neurons.z_neurons

