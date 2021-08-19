
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
% >> DataFiles = {'path/to/subj1.mat'; 'path/to/subj2.mat'};
% >> S = FNIRS1.DLM(DataFiles, 'McmcControl', fnirs1.mcmc_control);
%
% DataFiles specified should be load-able Matlab files that contain 
% variables that correspond to outcome data, named any of {'hbo', 'hbr',
% 'hbt'}; stimulus design, named 's'; and time stamps, named 't'.
%
% DataFiles can also be the empty character string,
% '', in which case FNIRS1.DLM will prompt the user to select participant
% data files. The call FNIRS1.DLM('') could be ambiguous with no other 
% arguments (does FNIRS1.DLM ask the user to select 'setup.dat' files? Or 
% participant data files?). In this case, FNIRS1.DLM will always resolve 
% this conflict by asking the user to select participant data files, and 
% passing these to fnirs1.specify_model with other defaults. If you need 
% help generating a list of 'setup.dat' files, FNIRS1 provides the 
% fnirs1.list_setup_files utility.
%
% Other valid Options are:  (see fnirs1.specify_model for full details)
%   'DownSampleRate'      - (numeric, integer)
%   'GroupData'           - (table)
%   'GroupFormula'        - (char)
%   'McmcControl'         - (fnirs1.mcmc_control)
%   'OutcomeType'         - (char)
%   'SpecificChannels'    - (numeric, integer)
%
%
% Example usage:
%   %% Single subject or group analysis, debug-length MCMC chains
%   fit = fnirs1.dlm('', 'DownSampleRate', 10, ...
%      'SpecificChannels', 1:2, ...
%      'McmcControl', fnirs1.mcmc_debug)
%
%   %% Group analysis: select multiple files at prompt,
%   % Existing demographic information contained in a table called 'Demo'
%   %   - Note: Names 'Cond' and 'TempDeriv' have special meanings and
%   %     should not appear in 'GroupData' table columns
%   fit = fnirs1.dlm('', 'DownSampleRate', 10, ...
%      'SpecificChannels', 1:2, ...
%      'GrouplData', Demo, 'GroupFormula', 'ID ~ Task * Cond', ...
%      'McmcControl', fnirs1.mcmc_control)
%   
%
% See also  fnirs1.specify_model, fnirs1.list_setup_files
%

summaries = [];
originaldir = pwd;

setups = parse_inputs(varargin{:});
success = false(size(setups));

% loop over setup files and fit models given in each
if (numel(setups) == 1)
    success = execute_setup(setups{1});
else
    parfor i = 1:numel(setups)
        success(i) = execute_setups(setups{i});
    end
end
cd(originaldir);  % execute_setups (below) changes directory

if (any(success))
    topsetupdir = fileparts(fileparts(setups{1}));
    if (contains(pwd, topsetupdir))
        cd(fullfile(topsetupdir, '..'));
        warning('Changing directory to %s', pwd);
    end
    
    % Next block:
    %  - loops over successful output files and read results summaries
    %  - zips directories with successful output files
    %  - removes un-zipped directories with successful output files
    loginfo = dir(fullfile(topsetupdir, '*', 'log', ...
        'Parameter_Estimates.log'));
    if ~isempty(loginfo)
        summaries = repmat(fnirs1.dlm_summary, size(loginfo));
        for i = 1:length(loginfo)
            % read summary info and set nicknames
            setupdir = fileparts(loginfo(i).folder);
            % summaries(i) = summaries(i).read_from_file(...
            %    fullfile(loginfo(i).folder, loginfo(i).name));
            % summaries(i).Nickname = fnirs1.utils.basename(setupDir);
            summaries(i) = summaries(i).parse_log_directory(...
                loginfo(i).folder);
            
            % zip successful folders and delete un-zipped versions
            try
                zip(setupdir, setupdir);
                % Remove the unzipped copy
                remove_folder_and_contents(setupdir);
            catch me
                warning(me.message);
            end
        end
    else
        warning('Cannot locate ''Parameter_Estimates.log'' file(s)');
        warning('Program may have terminated early. Check printed info');
    end
end
end



% Auxiliary functions
% -------------------

function setupfiles = parse_inputs(varargin)
% Parse input arguments. Forward to fnirs1.specify_model if necessary.
% Returns list of setup.dat files
if (nargin == 1)
    if isa(varargin{1}, 'fnirs1.modelspec')
        setupfiles = varargin{1}.setupfiles;
    elseif (ischar(varargin{1}) || isstring(varargin{1}))
        if (strcmp(varargin{1}, ''))
            % if empty char is the only input parameter, assume the user
            % wants to specify the model with specify_model('') (with
            % defaults)
            spec = fnirs1.specify_model('');
            setupfiles = spec.setupfiles;
        else
            setupfiles = cellstr(varargin{1});
        end
    elseif ~iscellstr(varargin{1})
        error("When called with with only one input, FNIRS1.DLM " + ...
            "takes either an fnirs1.modelspec argument or a " + ...
            "(set of) path(s) to 'setup.dat' file(s)");
    end
else
    spec = fnirs1.specify_model(varargin{:});
    setupfiles = spec.setupfiles;
end
if isempty(setupfiles)
    error('dlm:parse_inputs:emptysetup', ...
        'Cannot infer setup file location');
end
end


function success = fitdlm(setupfile)
% Create and execute SYSTEM fnrisdlm command
cmd = sprintf('%s %s', ...
    fullfile(fnirs1.home, 'include', 'fnirsdlm'), ...
    fnirs1.utils.basename(setupfile));
status = system(cmd, '-echo');
success = ~logical(status);
end


function success = execute_setup(setupfile)
if (exist(setupfile, 'file') ~= 2)
    error('''%s'' not found or is not a valid model setup file', ...
        setupfile);
end
setuppath = fileparts(setupfile);
cd(setuppath);
if (exist(fullfile(setuppath, 'log'), 'dir'))
    remove_folder_and_contents(fullfile(setuppath, 'log'));
    % ^^ fnirs1 private function
end
%
try success = fitdlm(setupfile);
catch me
    warning(me.identifier, 'fnirs1.dlm: caught error: %s', me.message);
end
end
