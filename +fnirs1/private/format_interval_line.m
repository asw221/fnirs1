
function interval = format_interval_line(line)
% Extract and format credible interval information from a character string.
% Return an empty vector if the line contains no relevant information, and
% return a struct with fields 'Level', 'Interval' otherwise
interval = [];
params = sscanf(line, '%f%% Cred.Int. = (%f, %f)');
if (~isempty(params))
    try
        interval = struct('Level', params(1) / 100, ...
            'Interval', params(2:3));
    catch ME
        warning('Unexpected formatting of credible interval line:\n\t%s', ...
            char(line));
    end
end
end
