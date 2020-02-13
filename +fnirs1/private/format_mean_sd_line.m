
function theta = format_mean_sd_line(line)
% Extract and format mean and SD information from a character string.
% Return empty vector if the line contains no relevant information, and
% return a struct with 'Estimate' and 'SE' fields otherwise
theta = [];
params = sscanf(line, 'mean = %f\tsd = %f');
if (~isempty(params))
    theta = struct('Estimate', params(1), 'SE', params(2));
end
end
