
function remove_folder_and_contents(folder)
if (~exist(folder, 'dir'))
    warning('Directory ''%s'' does not exist', folder);
else
    if (ispc)
        cmd = sprintf('rmdir %s /s /q', folder);
    else
        cmd = sprintf('rm -rf %s', folder);
    end
    system(cmd);
end
end
