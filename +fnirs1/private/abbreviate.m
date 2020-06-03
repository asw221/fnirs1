
function s = abbreviate(str)
% Abbreviate string-like data
s = regexprep(str, '[-_+~!@#$%^&*)(+=}{\[\]|\\><,.aeiou ]', '');
end
