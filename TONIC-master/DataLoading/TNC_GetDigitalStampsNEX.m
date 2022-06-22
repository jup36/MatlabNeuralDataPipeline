function [timestamps] = TNC_GetDigitalStampsNEX(datafile,id)
% FUNCTION DETAILS: Digital stamps need to be extracted somewhat differently from the nex and nev files. This is a utility function to assist with that.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

% index through the digital marker names (values.strings) and parse into
% separate variables looking for the specific id passed into the function.
counter =1;
disp(['Using id: ' id]);
if isfield(datafile,'markers')
    if isfield(datafile.markers{1}.values{1},'strings') % have to make sure this field exists. in some cases there are no events and the field is not written in nex file
        totalStrings = size(datafile.markers{1}.values{1}.strings,1);
        counter = 1;
        timestamps = 0;

        for i=1:totalStrings

            testForName = strcmp(datafile.markers{1}.values{1}.strings(i),id);
            if testForName==1
                indices(counter) = i;
                counter=counter+1;
            end

        end
    end
end

if counter>1
    % now run through the counter indices and retrieve and rescale the timestamps
    timestamps = datafile.markers{1}.timestamps(indices).*1000;
%     disp(sprintf('Total events (trials): %g',counter));
else
    timestamps = 0.001; 
%     disp('No events were found. Empty matrix returned.');
end
