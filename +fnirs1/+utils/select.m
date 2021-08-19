
function T = select(T, name, m)
% Select rows from table T by matching column NAME to input M
%
% Example usage:
%   lett = {'A'; 'B'; 'C'; 'D'; 'E'};
%   vals = randi(100, size(lett));
%   tab  = table(lett, vals);
%   fnirs1.utils.select(tab, 'lett', {'D', 'B', 'E'})
%   % Possible output:
%   % ans =
%   %   3×2 table
%   %     lett     vals
%   %     _____    ____
%   %     {'D'}     45 
%   %     {'B'}     19 
%   %     {'E'}     65 
%
if ~istable(T)
    error('First input to select must be a table');
end
%
% Verify that name is a valid variable name in T
name = char(name);
try T.(name);
catch ME
    rethrow(ME);
end
%
%
if iscellstr(m)
    m = string(m);  % == works for cellstr == string or vice versa
elseif iscell(m)
    error('Operator ''=='' is not supported for operands of type ''cell''');
end
%
% Execute subsetting
if islogical(m)
    T = T(m, :);  % Ignoring name
else
    % Find first matches for 'm' in column 'name'
    s = NaN(numel(m), 1);
    for i = 1:numel(m)
        f = find(T.(name) == m(i), 1, 'first');
        if ~isempty(f)
            s(i) = f;
        end
    end
    s = s(~isnan(s));
    T = T(s, :);
end
end
