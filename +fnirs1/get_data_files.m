
function files = get_data_files()
% FNIRS1.GET_DATA_FILES Open a dialogue window and request the user to 
% select file(s). Returns a cell array of full file paths
%
% Example usage:
%   files = fnirs1.get_data_files;
%
if ~(fnirs1.utils.interactive)
    error('%s cannot launch GUI for file selection', mfilename);
end
[files, path] = uigetfile(fullfile(last_path, '*'), ...
    'Select data file(s)', 'MultiSelect', 'on');
if (isnumeric(files))
    error("No data files selected");
end
last_path(path);  % save the path for later - user experience
files = cellstr(files)';
for i = 1:numel(files)
    files{i} = fullfile(path, files{i});
end
end
