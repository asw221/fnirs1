
function info = read_setup(fname)
% Read fnirs1 model setup file and return struct
subject_token = 'SUB_Replicates';
[names, values] = read_setup_components(fname);
[subinfo, substart] = parse_subject_info(names, values, subject_token);
if (substart > 0)
    values = [values(1:(substart-1)); {subinfo}];
    names = [names(1:(substart-1)); {'SubInfo'}];
end
info = cell2struct(values, names);
end

% ---

function [names, values] = read_setup_components(fname)
% Read the file, split raw text into name/value pairs
raw = fnirs1.utils.read_whole_file(fname);
values = cell(size(raw));
names = cell(size(raw));
keep = true(size(raw));
for i = 1:numel(raw)
    pair = strsplit(raw{i}, '=');
    if (numel(pair) == 2)
        names{i} = deblank(pair{1});
        v = strtrim(deblank(pair{2}));
        values{i} = str2double(v);
        if (isnan(values{i}) || isempty(values{i}))
            values{i} = v;
        end
    else
        keep(i) = false;
    end
end
values = values(keep);
names = names(keep);
end


function [c, i] = parse_subject_info(names, values, token)
% Parse subject-specific info into an nx1 cell array of cell arrays
starts = find(strcmpi(names, token));
n = numel(starts);
if (n > 0)
    starts = [starts; numel(names) + 1];
    c = cell(n, 1);
    for i = 1:n
        sub = starts(i):(starts(i+1) - 1);
        c{i} = [names(sub), values(sub)];
    end
    i = starts(1);
else
    c = {};
    i = -1;
end
end
