
function success = write_matrix(M, filename)
% Write matrix M to output file filename
success = false;
if (~isnumeric(M))
    error("Input M must be numeric-type");
end
if (~(ischar(filename) || isstring(output_path)))
    error("Input filename must be a valid path name");
end
fid = fopen(filename, 'w+');
if (fid ~= -1)
    [nrow, ncol] = size(M);
    for i = 1:nrow
        for j = 1:ncol
            if (j == ncol)
                fprintf(fid, '%f\n', M(i, j));
            else
                fprintf(fid, '%f ', M(i, j));
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
