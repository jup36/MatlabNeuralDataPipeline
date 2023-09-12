function datetimeObj = hmsToDateTime(hms)

% Define the time string
%hms = '14_58_14';

% Split the time string into hours, minutes, and seconds
timeParts = split(hms, '_');

% Convert the parts to numeric values
hours = str2double(timeParts{1});
minutes = str2double(timeParts{2});
seconds = str2double(timeParts{3});

% Create a datetime object with the specified time
datetimeObj = datetime('now', 'Format', 'HH:mm:ss'); % Get the current date with the desired time format
datetimeObj.Hour = hours;
datetimeObj.Minute = minutes;
datetimeObj.Second = seconds;

% You can now sort datetime objects
%datetimeArray = [datetimeObj]; % Add datetime objects to an array for sorting
%sortedDatetimeArray = sort(datetimeArray); % Sort the datetime array

end