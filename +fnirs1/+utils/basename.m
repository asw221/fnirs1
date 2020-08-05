
function filename = basename(fullFileName)
% Return the name of a file, removing its path
try
    [~, fname, fext] = fileparts(fullFileName);
    filename = sprintf('%s%s', fname, fext);
catch
    wasstring = isstring(fullFileName);
    try
        fullFileName = cellstr(fullFileName);
    catch ME
        error(ME.identifier, 'Input must be of a character type');
    end
    filename = cell(size(fullFileName));
    for i = 1:numel(fullFileName)
        [~, fname, fext] = fileparts(fullFileName{i});
        filename{i} = sprintf('%s%s', fname, fext);        
    end
    if wasstring
        filename = string(filename);
    end
end
end

