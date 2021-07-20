
function intervals = read_setup_interval_line(setupfile)
% Extract requested credible interval widths from setup file
intervals = [0.8, 0.9, 0.95, 0.99];
interval_token = 'Confidence_Intervals = ';
nit = numel(interval_token);
fid = fopen(setupfile, 'r');
if (fid ~= -1)
    token_found = false;
    while (~feof(fid) && ~token_found)
        tline = strtrim(fgetl(fid));
        if (numel(tline) >= nit && ...
                strcmp(tline(1:nit), interval_token))
            splt = strsplit(tline((nit+1):end), ' ');
            token_found = true;
        end
    end
    if token_found
        intervals = NaN(size(splt));
        ok = false(size(splt));
        for i = 1:numel(splt)
            intervals(i) = str2num(splt{i});
            ok(i) = isnumeric(intervals(i));
        end
        intervals = intervals(ok);
    end
else
    warning('Error opening %s', setupfile);
end
fclose(fid);
end
