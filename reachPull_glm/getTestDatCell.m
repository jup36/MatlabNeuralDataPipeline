function [testDatC] = getTestDatCell(dataSet,indices)
        for ii = 1:length(indices)
            testDatC{ii} = dataSet(:,indices(ii));
        end
end