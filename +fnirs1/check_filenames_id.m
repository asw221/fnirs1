
function [allok, varargout] = check_filenames_id(fnames, ids)
ok = contains(fnames, string(ids));
allok = all(ok);
if (nargout > 1)
    varargout{1} = find(~ok);
end
end
