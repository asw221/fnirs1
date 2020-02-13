
function setupFiles = list_setup_files(varargin)
% FNIRS1.LIST_SETUP_FILES List full file path instances of 'setup.dat'
% files in input directory, or any of the input directory's immediate
% children directories. Returns an empty cell array if no 'setup.dat' files
% are found
%
% FILES = FNIRS1.LIST_SETUP_FILES(DIR) returns a cell array of characters
% containing the paths of any 'setup.dat' files found. If DIR is not given,
% lists files for the current working directory is assumed
%
indir = pwd;
setupFiles = {};
if (nargin >= 1)
    indir = varargin{1};
end
if (nargin > 1)
    warning("FNIRS1.LIST_SETUP_FILES only takes one input argument. " + ...
        "Discarding the rest");
end
if (~exist(indir, 'dir'))
    error('%s not found', indir);
end

% Find files in indir
finfo = [ dir(fullfile(indir, 'setup.dat')); ...
    dir(fullfile(indir, '*', 'setup.dat')) ];
if (~isempty(finfo))
    setupFiles = cell(length(finfo), 1);
    for i = 1:length(finfo)
        setupFiles{i} = fullfile(finfo(i).folder, finfo(i).name);
    end
end
end

