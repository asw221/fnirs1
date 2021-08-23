
function success = write_matrix(M, filename, varargin)
% Write matrix M to output file filename
% Optional third argument specifies the number of decimal places to write.
% Default is 15
success = false;
if (~isnumeric(M))
    error("Input M must be numeric-type");
end
if (~(ischar(filename) || isstring(filename)))
    error("Input filename must be a valid path name");
end
if (nargin <= 2)
    digits = 15;
else
    digits = max(fix(varargin{1}), 1);
end
fid = fopen(filename, 'w+');
fmt = sprintf('%%.%if', digits);
if (fid ~= -1)
    [nrow, ncol] = size(M);
    for i = 1:nrow
        for j = 1:ncol
            fprintf(fid, fmt, M(i, j));
            if (j == ncol)
                fprintf(fid, '\n');
            else
                fprintf(fid, ' ');
            end
        end
        % fprintf(fid, '%s\n', sprintf('%f ', M(i, :)));
    end
    success = true;
else
    warning("Error opening " + string(filename));
end
fclose(fid);
end
