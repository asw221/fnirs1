
function filename = basename(fullFileName)
% return the name of a file, removing the path to it
[~, fname, fext] = fileparts(fullFileName);
filename = sprintf('%s%s', fname, fext);
end

