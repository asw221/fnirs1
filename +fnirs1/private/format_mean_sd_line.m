
function theta = format_mean_sd_line(line)
% Extract and format mean and SD information from a character string.
% Return empty vector if the line contains no relevant information, and
% return a struct with 'Estimate' and 'SE' fields otherwise
theta = struct();
% params = sscanf(line, 'mean = %f\tsd = %f');
name = deblank(strsplit(line, 'mean = '));
if (numel(name) == 1)
    name = {''};
end
params = regexp(line, '\d+\.?\d*', 'match');
for i = 1:numel(params)
    params{i} = sscanf(params{i}, '%f');
end
params = [params{:}];
N = numel(params);
if (~isempty(params))
    theta = struct('Name', name{1}, 'Estimate', params(N-1), ...
        'SE', params(N));
end
end
