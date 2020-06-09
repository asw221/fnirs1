
function s = abbreviate(str)
% Abbreviate string-like data. Remove lower case vowels and some symbol
% characters
s = regexprep(str, '[~!@#$%^&)(}{\[\]|\\,aeiou]', '');
end
