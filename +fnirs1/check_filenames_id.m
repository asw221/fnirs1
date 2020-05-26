
function [allok, varargout] = check_filenames_id(fnames, ids)
% FNIRS1.CHECK_FILENAMES_ID  Perform a quick check to make sure participant
% ID's appear in a desired order
%
% This check is NOT foolproof, but may help to catch and correct mistakes
% when dealing with a large number of files.
%
% ok = FNIRS1.CHECK_FILENAMES_ID(fnames, ids) returns logical true if there
% are simple substrings in fnames that match ids, and they appear in the
% same order. This may be a useful check if data file names contain
% participant identifiers.
%
% [ok, mismatch] = FNIRS1.CHECK_FILENAMES_ID(fnames, ids) additionally
% returns an integer index (mismatch) of elements in fnames that
% do not match the order of ids. It may be important to note that fnames is
% tested specifically against ids: if all N fnames match the first N ids,
% then FNIRS1.CHECK_FILENAMES_ID will return true and an empty mismatch.
%
% Example usage:
%   data_files = {'path/to/subj1.mat'; 'path/to/subj2.mat'}
%   ok = fnirs1.check_filenames_id(data_files, [1, 2])
%   [ok2, msmtch] = fnirs1.check_filenames_id(data_files, {'0', '2'})
%
if ~isstring(fnames)
    try
        fnames = string(fnames);
    catch ME
        error(ME.message);
    end
end
if ~isstring(ids)
    try
        ids = string(ids);
    catch ME
        error(ME.message);
    end
end
ok = contains(fnames, ids);
allok = all(ok);
if (nargout > 1)
    varargout{1} = find(~ok);
end
end
