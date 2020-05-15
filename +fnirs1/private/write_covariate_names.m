
function success = write_covariate_names(nms, varargin)
% WRITE_COVARIATE_NAMES  Write a 'names.txt' file to the input directory.
%   S = WRITE_COVARIATE_NAMES(nms) writes nms to a file in the current 
%   directory. Input nms should be a cellstr or array of strings. 
%   Returns boolean S (true if successful)
%
%   S = WRITE_COVARIATE_NAMES(nms, dirname) Writes nms to a file in 
%   directory dirname and returns boolean S (true if successful)
%
success = false;
if ~(iscellstr(nms) || isstring(nms))
    error("write_covariate_names requires cellstr or string array input");
end
nms = cellstr(nms);
output_path = pwd;
if (nargin > 1)
    output_path = varargin{1};
    if (nargin > 2)
        warning("discarding " + (nargin - 1) + " unused arguments");
    end
end
if (~(ischar(output_path) || isstring(output_path)))
    error("Input must be a valid path name");
end
output_file = fullfile(output_path, 'names.txt');
fid = fopen(output_file, 'w+');
if (fid ~= -1)
    for i = 1:numel(nms)
        fprintf(fid, '%s\n', nms{i});
    end
    success = true;
else
    warning("Error opening " + convertCharsToStrings(output_file));
end
fclose(fid);
end
