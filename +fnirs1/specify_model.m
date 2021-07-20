
function spec = specify_model(datafiles, varargin)
% FNIRS1.SPECIFY_MODEL Reads input data and parameters and writes a
% temporary directory containing all the files needed to run
% subject-specific or group-level analyses. Returns an fnirs1.modelspec
% object
%
% DIR = FNIRS1.SPECIFY_MODEL('') opens a graphical interface to select data
% files; other parameters set to default values
%
% DIR = FNIRS1.SPECIFY_MODEL(FILES) writes the data in FILES to DIR; other
% parameters set to default values
%
% DIR = FNIRS1.SPECIFY_MODEL(FILES, ParamName, ParamValue, ...) sets
% specific parameters. See below for valid options
%
% In general, FILES can either be a empty, numeric, the path to a 
% single file, or a cell array of paths to files. If empty or numeric,
% FNIRS1.SPECIFY_MODEL will open a system interface allowing the user to
% select files manually. In case of numeric input, the argument will be
% interpreted to select per participant (useful for repeated measures on
% the same subject). In this case, when the file selection GUI opens, users
% should select files for all patients for replicate 1 and close the GUI.
% Subsequent GUIs will open for the remaining replicates. Users must take
% care to select the same number of files for each replicate: fnirs1 does
% not support imbalanced designs.
% GUI-based file selection is of course only available when running Matlab
% in interactive mode.
%
% If FILES is given as a single file, the output data in
% DIR will be written to specify a single-subject analysis. If multiple 
% files are given, the setup files will be written to specify a group-level
% analysis. For group-level analyses with repeated measures, FILES should
% be provided as a rectangular cell array or matrix of strings, where each
% row i corresponds to set of files for patient i, and each column j lists
% the files for the j-th repeated measure.
% In the above cases, analyses are set up to be done independently
% over all channels; subsets of channels can be specified (see
% 'SpecificChannels' option below).
%
% FILES specified should be load-able Matlab files that contain 
% variables that correspond to outcome data, named any of {'hbo', 'hbr',
% 'hbt'}; stimulus design, named 's'; and time stamps, named 't'.
%
% Optional parameters can be set using ParamName, ParamValue pairs. Valid
% options are:
%   'CredibleIntervals' - (numeric) a vector of credible interval widths
%       to output. Default is the vector [0.8, 0.9, 0.95, 0.99] for 80%,
%       ..., 99% posterior credible intervals
%
%   'DownSampleRate' - (numeric) an integer factor of participant sampling
%       rates. Should be length 1 or a vector the same length as the number
%       of data FILES. Default is 1 which does not perform down-sampling
%
%   'GroupData' - (table) demographic information for group
%       analyses. There are two special/reserved variable names that should
%       not appear in a 'GroupData' table directly. They are: 'Cond' for
%       task condition (derived from the columns of participants' task
%       design matrix), and 'TempDeriv' for HRF temporal derivatives (if
%       requested). FNIRS1 will always append 'Cond' as a column to group
%       data; 'TempDeriv' is optionally added as a main effect directly to
%       the group level covariate matrix. This argument is required for
%       group-level analyses
%
%   'GroupFormula' - (char) will look something like, 'y ~ age + sex'.
%       While the formula requires a numeric outcome variable (like 'y' 
%       above), this can be any continuous variable in 'GroupData' and will
%       not actually be referred to in model fitting. (See 'GroupData'
%       above) The special variable 'Cond' can be referred to in the group
%       formula; 'TempDeriv' should not be. If using 'Cond' in a model 
%       formula, please make sure it appears LAST of all the categorical
%       covariates referred to. For example, 'y ~ Task * Cond' is ok, but
%       'y ~ Cond + Task' is not (assuming Task is categorical).
%       This argument is required for group-level analyses
%
%   'McmcControl' - (fnirs1.mcmc_control) can be any valid
%       fnirs1.mcmc_control object. Default is the same as returned by
%       fnris1.mcmc_control
%
%   'OutcomeType' - (char) should be one of {'hbo', 'hbr', 'hbt'}, to
%       specify which data to use for the outcome variable. Default is
%       'hbo'
%
%   'SpecificChannels' - (numeric, integer) should be an integer vector to
%       index specific channels if it's desired to perform an analysis over
%       subsets of channels. Default includes all channels in the analyses
%       and does not subset
%
% Example usage:
%   %% Single subject or group analysis, debug-length MCMC chains
%   setup_dir = fnirs1.specify_model('', 'DownSampleRate', 10, ...
%      'SpecificChannels', 1:2, ...
%      'McmcControl', fnirs1.mcmc_debug)
%   setup_files = fnirs1.list_setup_files(setup_dir)
%   fit = fnirs1.dlm(setup_files)   
%
%   %% Group analysis: select multiple files at prompt,
%   % Existing demographic information contained in a table called 'Demo'
%   %   - Note: Names 'Cond' and 'TempDeriv' have special meanings and
%   %     should not appear in 'GroupData' table columns
%   setup_dir = fnirs1.specify_model('', 'DownSampleRate', 10, ...
%      'GrouplData', Demo, 'GroupFormula', 'ID ~ Task * Cond', ...
%      'McmcControl', fnirs1.mcmc_control)
% 
%
% See also
% fnirs1.dlm, fnirs1.modelspec, fnirs1.list_setup_files,
% fnirs1.mcmc_control, fnirs1.mcmc_debug, load, readtable
% 

if (isempty(datafiles) || isnumeric(datafiles))
    datafiles = select_data_files(datafiles);
else
    try datafiles = cellstr(datafiles);
    catch me
        error('specify_model:improperinput', ...
            ['Argument ''datafiles'' should be a scalar or a type ', ...
            'convertable to cellstr']);
    end
end


spec = fnirs1.modelspec(datafiles, varargin{:});

fprintf('Reading patient data and writing temporaries\n');
fprintf('This may take several minutes...\n');

spec = spec.write_channel_data();

fprintf('\tModel setup complete.\n\n');
end



% Auxiliary functions
% -------------------

function datafiles = select_data_files(code)
% Open fnirs1 file selection GUI if able and prompt the user to select data
% files. Returns a cell array of strings
if ~fnirs1.utils.interactive
    error('specify_model:notinteractive', ...
        'Interactive file selection not available');
end
%
nrep = 1;
if isnumeric(code)
    nrep = max(fix(code), 1);
end
%
df = fnirs1.get_data_files();
datafiles = cell(numel(df), 1);
datafiles(:, 1) = df;
if (nrep > 1)
    for r = 2:nrep
        df = fnirs1.get_data_files();
        if (numel(df) ~= size(datafiles, 1))
            error('specify_model:select_data_files:imbalanceddesign', ...
                ['Inconsistent numbers of files selected. fnirs1 ', ...
                'does not support imbalanced designs']);
        end
        datafiles(:, r) = df;
    end
end
end
