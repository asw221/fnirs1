
function TF = regexpl(str, pattern)
% Logical matches for regular expressions
str = cellstr(str);
pattern = char(pattern);
TF = false(size(str));
for i = 1:numel(str)
    val = regexp(str{i}, pattern, 'ONCE');
    TF(i) = ~isempty(val);
end
end
