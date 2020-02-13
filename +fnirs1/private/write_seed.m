
function success = write_seed(varargin)
% WRITE_SEED  Write a 'seed.dat' file to the input directory.
%   S = WRITE_SEED() writes a seed file to the current directory and
%   returns boolean S (true if successful)
%
%   S = WRITE_SEED(DIRNAME) Writes a 'seed.dat' file to DIRNAME and
%   returns boolean S (true if successful)
%
% See also RNG
%
success = false;
if (nargin > 0)
    output_path = varargin{1};
    if (nargin > 1)
        warning("discarding " + (nargin - 1) + " unused arguments");
    end
else
    output_path = pwd;
end
if (~(ischar(output_path) || isstring(output_path)))
    error("Input must be a valid path name");
end
output_file = fullfile(output_path, 'seed.dat');
fid = fopen(output_file, 'w+');
if (fid ~= -1)
    fprintf(fid, '%i%i %i%i %i', randi([0, intmax], 5, 1));
    success = true;
else
    warning("Error opening " + convertCharsToStrings(output_file));
end
fclose(fid);
end
