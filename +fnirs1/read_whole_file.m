
function lines = read_whole_file(fname)
% Reads whole file into cell array: one cell for each line
lines = {};
fid = fopen(fname, 'r');
if (fid ~= -1)
    tline = strtrim(fgetl(fid));
    lines = {tline};
    while ~feof(fid)
        tline = strtrim(fgetl(fid));
        lines = [lines; tline];
    end
else
    warning('Error opening %s', fname);
end
fclose(fid);
end
