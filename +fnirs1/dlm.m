
function summaries = dlm(varargin)
% FNIRS1.DLM Fit Dr. Johnson's distributed lag-type models to fNIRs data
%
% S = FNIRS1.DLM(INFILE) Fits the model specified by INFILE, which should
% be of a particular layout (see example setup.dat's). Returns an
% fnirs1.dlm_summary object (or vector of objects). INFILE can also be 
% given as a character cell array with paths to different setup files in 
% each cell.
%
% S = FNIRS1.DLM(DataFiles, Option, OptionValue, ...) may be more
% convenient for first use. In this case, DataFiles and all optional
% arguments will be forwarded to fnirs1.specify_model, which will run
% before model fitting begins. Note that the user must specify at least one
% Option to invoke this usage, as when called with one argument, FNIRS1.DLM
% assumes the setup.dat-based syntax above (with one exception noted 
% below). One easy way to specify an optional argument without changing any
% defaults is to call, for example,
% >> S = FNIRS1.DLM(DataFiles, 'McmcControl', fnirs1.mcmc_control);
%
% Note that the call FNIRS1.DLM('') is ambiguous with no other arguments 
% (does FNIRS1.DLM ask the user to select 'setup.dat' files? Or participant
% data files?). In this case, FNIRS1.DLM will always resolve this conflict 
% by asking the user to select participant data files, and passing these to
% fnirs1.specify_model with other defaults. If you need help generating a
% list of 'setup.dat' files, FNIRS1 provides the fnirs1.list_setup_files
% utility.
%
% See also  fnirs1.specify_model, fnirs1.list_setup_files
%

summaries = [];

% Organize list of setup files or call specify_model to do this
% automatically
if (nargin == 1)
    setupFiles = varargin{1};
    if (ischar(setupFiles) || isstring(setupFiles))
        if (strcmp(setupFiles, ''))
            % if empty char is the only input parameter, assume the user
            % wants to specify the model with specify_model('') (with
            % defaults)
            setupFiles = fnirs1.list_setup_files(fnirs1.specify_model(''));
        else
            setupFiles = cellstr(setupFiles);
        end
    elseif ~iscellstr(setupFiles)
        error("When called with with only one input, FNIRS1.DLM " + ...
            " assumes the argument is a (set of) path(s) to " + ...
            "'setup.dat' file(s)");
    end
else
    setupFiles = fnirs1.list_setup_files(fnirs1.specify_model(varargin{:}));
end

success = false(size(setupFiles));

% loop over setup files and fit models given in each
original_dir = pwd;
for i = 1:length(setupFiles)
    if (exist(setupFiles{i}, 'file') ~= 2)
        error('''%s'' not found or is not a valid model setup file', ...
            setupFiles{i});
    end
    setupPath = fileparts(setupFiles{i});
    cd(setupPath);
    if (exist(fullfile(setupPath, 'log'), 'dir'))
        remove_folder_and_contents(fullfile(setupPath, 'log'));
    end
    try 
        fnirs1.fitDlm(basename(setupFiles{i}));
        success(i) = true;
    catch ME
        warning(ME.identifier, 'fnirs1.dlm: caught error: %s', ME.message);
    end
end

cd(original_dir);
if (any(success))
    topSetupDir = fileparts(fileparts(setupFiles{1}));
    if (contains(pwd, topSetupDir))
        cd(fullfile(topSetupDir, '..'));
        warning('Changing directory to %s', pwd);
    end
    
    % Next block:
    %  - loops over successful output files and read results summaries
    %  - zips directories with successful output files
    %  - removes un-zipped directories with successful output files
    loginfo = dir(fullfile(topSetupDir, '*', 'log', 'Parameter_Estimates.log'));
    if ~isempty(loginfo)
        summaries = repmat(fnirs1.dlm_summary, size(loginfo));
        for i = 1:length(loginfo)
            % read summary info and set nicknames
            setupDir = fileparts(loginfo(i).folder);
            summaries(i) = summaries(i).read_from_file(...
                fullfile(loginfo(i).folder, loginfo(i).name));
            summaries(i).Nickname = basename(setupDir);
            
            % zip successful folders and delete un-zipped versions
            zip(setupDir, setupDir);
            remove_folder_and_contents(setupDir);
        end
    else
        warning('Cannot locate ''Parameter_Estimates.log'' file(s)');
    end
end
end
