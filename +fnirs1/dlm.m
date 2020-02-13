
function summaries = dlm(varargin)
% DLM Fit Dr. Johnson's distributed lag-type model to fNIRs data.
%
% S = FNIRS1.DLM(INFILE) Fits the model specified by INFILE, which should
% be of a particular layout (see documentation). Returns boolean S (true if
% model was fit successfully)
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
        error('%s not found', setupFiles{i});
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
