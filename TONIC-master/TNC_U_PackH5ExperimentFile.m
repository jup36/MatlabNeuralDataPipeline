function fileName = TNC_U_PackH5ExperimentFile(blackrock_file_path, behavior_data_path, preprocess_flag, nev_flag)

Write the data to the HDF5 file.
h5write('my_example_file.h5', '/dataset1', testdata)

%%
if nev_flag

    h5create(fileName, '/blackrock/nev', size(nev_data));
    
end
if ns2_true

        h5create(fileName, '/blackrock/ns2', size(ns2_data));

end
if ns3_true

        h5create(fileName, '/blackrock/ns3', size(ns3_data));

end
if ns4_true

        h5create(fileName, '/blackrock/ns4', size(ns4_data));

end

if ns5_true

        h5create(fileName, '/blackrock/ns5', size(ns5_data));

end



%%
if numel(behavior_data_path)>0

        h5create(fileName, '/behavior', size(behav_data));


end

%%
if preprocess_flag
    
    
end