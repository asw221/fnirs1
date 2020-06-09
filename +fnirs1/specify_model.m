
function outdir = specify_model(dataFiles, varargin)
% FNIRS1.SPECIFY_MODEL Reads input data and parameters and writes a
% temporary directory containing all the files needed to run
% subject-specific or group-level analyses. Returns the name of that
% temporary directory
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
% In general, FILES can either be an empty character, '', the path to a 
% single file, or a cell array of paths to files. If empty,
% FNIRS1.SPECIFY_MODEL will open a system interface allowing the user to
% select files manually. If only a single file is given, the output data in
% DIR will be written to specify a single-subject analysis. If multiple 
% files are given, the setup files will be written to specify a group-level
% analysis. In both cases, analyses are set up to be done independently
% over all channels; subsets of channels can be specified.
%
% FILES specified should be load-able Matlab files that contain 
% variables that correspond to outcome data, named any of {'hbo', 'hbr',
% 'hbt'}; stimulus design, named 's'; and time stamps, named 't'.
%
% Optional parameters can be set using ParamName, ParamValue pairs. Valid
% options are:
%   'DownSampleRate' - (numeric) an integer factor of participant sampling
%       rates. Should be length 1 or a vector the same length as the number
%       of data FILES. Default is 1 which does not perform down-sampling
%
%   'GroupCovariates' - (numeric; optional) a covariates matrix for group
%       analysis. Should have the same number of rows as the number of data
%       FILES. Default is empty for single-subject analyses, or a vector of
%       1's for group-level analyses. Instead of using this parameter/value
%       pair directly, however, we recommend that users should typically 
%       specify the 'GroupData' and 'GroupFormula' options instead
%
%   'GroupCovariateNames' - (cellstr; optional) to be used mainly in
%       conjunction with the 'GroupCovariates' option. These can also be
%       inferred from syntax using 'GroupData' and 'GroupFormula' 
%
%   'GroupData' - (table; optional) demographic information for group
%       analyses. There are two special/reserved variable names that should
%       not appear in a 'GroupData' table directly. They are: 'Cond' for
%       task condition (derived from the columns of participants' task
%       design matrix), and 'TempDeriv' for HRF temporal derivatives (if
%       requested). FNIRS1 will always append 'Cond' as a column to group
%       data; 'TempDeriv' is optionally added as a main effect directly to
%       the group level covariate matrix
%
%   'GroupFormula' - (char; optional unless 'GroupData' is provided) this
%       will look something like, 'y ~ age + sex'. While the formula
%       requires a numeric outcome variable (like 'y' above), this can be
%       any continuous variable in 'GroupData' and will not actually be
%       referred to in model fitting. (See 'GroupData' above) The special
%       variable 'Cond' can be referred to in the group formula;
%       'TempDeriv' should not be. If using 'Cond' in a model formula,
%       please make sure it appears LAST of all the categorical
%       covariates referred to. For example, 'y ~ Task * Cond' is ok, but
%       'y ~ Cond + Task' is not (assuming Task is categorical).
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
% fnirs1.dlm, fnirs1.list_setup_files, fnirs1.mcmc_control, 
% fnirs1.mcmc_debug, load, readtable
% 


options = struct(...
    'DownSampleRate', 1, ...
    'GroupCovariates', [], ...
    'GroupCovariateNames', {''}, ...
    'GroupData', table(), ...
    'GroupFormula', '', ...
    'McmcControl', fnirs1.mcmc_control, ...
    'OutcomeType', 'hbo', ...
    'SpecificChannels', [] ...
    );
if (fix(length(varargin) / 2) ~= length(varargin) / 2)
    error('FNIRS1.SPECIFY_MODEL uses property name/value pairs');
end

% Parse varargin
for pair = reshape(varargin, 2, [])
    nm = pair{1};
    if (any(strcmp(nm, fieldnames(options))))
        checkParameterPair(pair);
        options.(nm) = pair{2};
    else
        error('%s is not a recognized parameter', nm);
    end
end
if (~ischar(options.GroupFormula))
    options.GroupFormula = char(options.GroupFormula);
end
if (~ischar(options.OutcomeType))
    options.OutcomeType = char(options.OutcomeType);
end
if (~isempty(options.SpecificChannels))
    options.SpecificChannels = sort(unique(floor(options.SpecificChannels)));
end



% Check dataFiles argument
if (isnumeric(dataFiles) || isempty(dataFiles))
    dataFiles = fnirs1.get_data_files;
end
if (ischar(dataFiles) || isstring(dataFiles))
    dataFiles = cellstr(dataFiles);
elseif (~iscellstr(dataFiles))
    error("dataFiles argument must be UI flag (empty) " + ...
        "or convertable to a cell-array of strings");
end


% Adjust set parameters
N = length(dataFiles);
groupAnalysis = N > 1;
data = load(dataFiles{1});
M = size(data.s, 2);        % Number of task conditions


%  - options.DownSampleRate
if (numel(options.DownSampleRate) == 1)
    options.DownSampleRate = repmat(options.DownSampleRate, N, 1);
elseif (numel(options.DownSampleRate) ~= N)
    error("Length of participant down-sampling rates should either be 1 " + ...
        "or the same length as dataFiles");
end

%  - options.GroupCovariates
if (groupAnalysis)
    % Process GroupFormula and GroupData if given. GroupData requires a
    % GroupFormula
    if (isempty(options.GroupFormula) && ~isempty(options.GroupData))
        error('''GroupData'' parameter must be accompanied by a ''GroupFormula'' parameter');
    elseif (~isempty(options.GroupFormula) && isempty(options.GroupData))
        warning('''GroupFormula'' will be ignored without ''GroupData''');
    end
    if (~isempty(options.GroupData))
        % remove any instance of TempDeriv from the formula
        if (contains(options.GroupFormula, 'TempDeriv'))
            % regular expression is not foolproof (tertiary+ interactions)
            fprintf('%s -> ', options.GroupFormula);
            options.GroupFormula = regexprep(options.GroupFormula, ...
                '[ *+-(:]*TempDeriv[):]*', '');
            fprintf('%s\n', options.GroupFormula);
            warning('Don''t put TempDeriv directly in model formula');
        end
        % Add Cond column to GroupData
        options.GroupData = fnirs1.expand_table_conditions(...
            options.GroupData, M, 'Cond');
        % Parse group covariates matrix
        design = fnirs1.mixed_effects_design(options.GroupData, ...
            options.GroupFormula);
        design = design.reformat_base_condition;
        options.GroupCovariates = designMatrix(design, 'Fixed');
        options.GroupCovariateNames = design.CoefficientNames;
    end
end
if (isempty(options.GroupCovariates))
    options.GroupCovariateNames = {'(Intercept)'};
    options.GroupCovariates = ones(N * M, 1);
    % options.GroupCovariates = eye(N);
elseif (size(options.GroupCovariates, 1) ~= N * M)
    error("'GroupCovariates' must have number of rows equal to " + ...
        "the number of participants/data files times the " + ...
        "number of columns in the task design");
end

% Add Temporal derivatives if requested
if (options.McmcControl.includeDerivatives)
    % Inject TempDeriv columns into GroupCovariates and append
    % GroupCovariateNames
    nx = size(options.GroupCovariates, 1);
    P = fix(nx / N);  % Should be = M 99% of the time
    
    ndx = [reshape(1:nx, P, N)', repmat((1:M) + nx, N, 1)]';
    ndx = reshape(ndx, numel(ndx), 1);
    options.GroupCovariates = blkdiag(...
        options.GroupCovariates, eye(M));
    options.GroupCovariates = options.GroupCovariates(ndx, :);
    options.GroupCovariateNames = [options.GroupCovariateNames, ...
        cellstr("TempDeriv_" + string(1:M))];
end

% Setup output directory and write data 
outdir = fullfile(pwd, sprintf('_fnirs_%s_data_', ...
    datetime('now', 'TimeZone', 'local', 'Format', 'd-MMM-y-HH-mm-ss-SSS')));
if (~exist(outdir, 'dir'))
    mkdir(outdir);
end

% make sub-directories for each channel and load -> write data
nChannels = 0;
setupFiles = {};
fprintf("Reading participant data files and writing temporaries.\n" + ...
    "This may take a minute\n");
for i = 1:N
    % check that each data file exists and has the correct fields
    [~, participantFileBaseName] = fileparts(dataFiles{i});
    try
        data = load(dataFiles{i});
        if (size(data.s, 2) ~= M)
            error('Mismatch in number of task conditions: ''%s'' and ''%s''', ...
                dataFiles{1}, dataFiles{i});
        end
    catch ME
        remove_folder_and_contents(outdir);
        error(ME.identifier, '%s', ME.message);
    end
    msng = findMissingFields(data, options.OutcomeType);
    if (~isempty(msng))
        remove_folder_and_contents(outdir);
        error('File ''%s''\n\tis missing fields:  %s\b\b', ...
            dataFiles{i}, sprintf('%s, ', string(msng)));
    end
    
    outcome = data.(options.OutcomeType);
    samplingRate = floor(1 / (data.t(2) - data.t(1)));
    if (~validDownSampleRate(samplingRate, options.DownSampleRate(i)))
        remove_folder_and_contents(outdir);
        error('Improper down-sampling rate (%.2f) for data from file:\n\t%s', ...
            options.DownSamplerate(i), dataFiles{i});
    end
    
    if (size(outcome, 1) ~= size(data.s, 1))
        remove_folder_and_contents(outdir);
        error('File %s timepoint mismatch:\n\tField ''%s'' has %d entries while ''s'' has %d', ...
            fnirs1.utils.basename(dataFiles{i}), options.OutComeType, size(outcome, 1), ...
            size(data.s, 1));
    end
    
    % If we ever need to append S matrix for TempDeriv, that will happen
    % here, for example, like this:
    % % If TempDeriv option requested, append I_M block to task design matrix
    % if (options.McmcControl.includeDerivatives)
    %     data.s = blkdiag(data.s, eye(M));
    %     outcome = [outcome; zeros(M, size(outcome, 2))];  % <- 
    % end
    
    
    if (i == 1)
        nChannels = size(outcome, 2);
        setupFiles = cell(nChannels, 1);
        if (isempty(options.SpecificChannels))
            options.SpecificChannels = 1:nChannels;
        elseif (any(options.SpecificChannels > nChannels))
            remove_folder_and_contents(outdir);
            error('Input ''SpecificChannels'' outside range of channels in data');
        end
        
        for ch = options.SpecificChannels
            chdir = fullfile(outdir, sprintf('ch%04d', ch));
            if (~exist(chdir, 'dir'))
                mkdir(chdir);
            end
            setupFiles{ch} = fullfile(chdir, 'setup.dat');
            
            % Create setup.dat's and write their preambles
            setupFid = fopen(setupFiles{ch}, 'w');
            if (setupFid == -1)
                remove_folder_and_contents(outdir);
                error('Cound not write %s', setupFiles{ch});
            end
            writeSetupPreamble(setupFid, N, data.s, options.McmcControl);
            fclose(setupFid);
            
            % Write seed.dat, names.txt & covar.dat
            if (~write_seed(chdir))
                remove_file_and_contents(outdir);
                error('Could not write seed.dat');
            end
            if (~write_covariate_names(options.GroupCovariateNames, ...
                    chdir))
                remove_file_and_contents(outdir);
                error('Could not write names.txt');
            end
            if (~write_matrix(options.GroupCovariates, ...
                    fullfile(chdir, 'covar.dat')))
                remove_file_and_contents(outdir);
                error('Could not write covar.dat');
            end
            
        end  % for ch = options.SpecificChannels
    end  % if (i == 1)
    
    if (size(outcome, 2) ~= nChannels)
        remove_folder_and_contents(outdir);
        error('Number of channels mismatch between participants:\n\t%s', ...
            dataFiles{i});
    end
    
    % Loop back over channels and write participant-specific data
    for ch = options.SpecificChannels
        chdir = fileparts(setupFiles{ch});
        participantDataFile = fullfile(chdir, ...
            sprintf('%s_%s.txt', participantFileBaseName, options.OutcomeType));
        participantDesignFile = fullfile(chdir, ...
            sprintf('%s_stim.txt', participantFileBaseName));
        
        % write out data
        if (~write_matrix(outcome(:, ch), participantDataFile))
            remove_folder_and_contents(outdir);
            error('Could not write %s', participantDataFile);
        end
        if (~write_matrix(data.s, participantDesignFile))
            remove_folder_and_contents(outdir);
            error('Could not write %s', participantDesignFile);
        end
        
        % append subject info to end of channel setup file
        %  - at this point the file has already been opened and written to
        %    earlier, so should be successful
        setupFid = fopen(setupFiles{ch}, 'a');
        writeSetupParticipant(setupFid, participantDataFile, ...
            participantDesignFile, samplingRate, ...
            options.DownSampleRate(i));
        fclose(setupFid);
    end  % for ch = options.SpecificChannels
    fprintf('\tFiles written for: %s\n', participantFileBaseName);
end  % for i = 1:N
fprintf('Model setup complete\n');
end






% Auxiliary Functions
% --------------------------------------------------------------------


function checkParameterPair(pair)
% Ensure input parameters in pair (cell) contain appropriate types    
% 'GroupCovariateNames', {''}, ...
%     'GroupData', table(), ...
%     'GroupFormula', '', ...
nm = pair{1};
if (strcmp(nm, 'DownSampleRate') && ~isnumeric(pair{2}))
    error('''DownSampleRate'' should be a numeric parameter');
end
if (strcmp(nm, 'GroupCovariates') && ~isnumeric(pair{2}))
    error('''GroupCovariates'' should be a numeric parameter');
end
if (strcmp(nm, 'GroupCovariateNames') && ...
        ~(iscellstr(pair{2}) || isstring(pair{2})))
    error('''GroupCovariateNames'' should be a cell array of strings');
end
if (strcmp(nm, 'GroupData') && ~istable(pair{2}))
    error('''GroupData'' should be a table object');
end
if (strcmp(nm, 'GroupFormula') && ~(ischar(pair{2}) || isstring(pair{2})))
    error('''GroupFormula'' should be a character string');
end
if (strcmp(nm, 'McmcControl') && ~isa(pair{2}, 'fnirs1.mcmc_control'))
    error('''McmcControl'' should be an fnirs1.mcmc_control object');
end
if (strcmp(nm, 'OutcomeType'))
    if ~(ischar(pair{2}) || isstring(pair{2}))
        error('''OutcomeType'' should be a char/string parameter');
    end
    if ~any(strcmpi(pair{2}, {'hbo', 'hbr', 'hbt'}))
        error('''OutcomeType'' should be one of {''hbo'', ''hbr'', ''hbt''}');
    end
end
if (strcmp(nm, 'SpecificChannels'))
    if (~isnumeric(pair{2}))
        error('''SpecificChannels'' should be a numeric parameter');
    elseif (any(pair{2} < 1))
        error('''SpecificChannels'' should all be positive integers');
    end
end
end




function missingFields = findMissingFields(data, outcomeType)
% Return names of any missing fields in data
fn = fieldnames(data);
necessaryFields = {outcomeType; 's'; 't'};
ok = logical(size(necessaryFields));
for i = 1:numel(necessaryFields)
    ok(i) = any(strcmp(necessaryFields{i}, fn));
end
missingFields = necessaryFields(~ok);
end




function success = validDownSampleRate(samplingRate, downSampleRate)
% Return true if downSampleRate is a positive integer factor of
% samplingRate
success = (downSampleRate > 0) && ...
    (fix(samplingRate / downSampleRate) == samplingRate / downSampleRate);
end




function success = writeSetupParticipant(fid, dataFile, designFile, ...
    samplingRate, downSampleRate)
% Write subject-specific information to fid
success = true;
fprintf(fid, '\nSUB_Replicates = %d\n', 1);
fprintf(fid, 'SUB_Data = %s\n', fnirs1.utils.basename(dataFile));
fprintf(fid, 'SUB_Design = %s\n', fnirs1.utils.basename(designFile));
fprintf(fid, 'SUB_Freq = %d\n', samplingRate);
fprintf(fid, 'SubSamp_Freq = %d\n', downSampleRate);
end




function success = writeSetupPreamble(fid, N, design, mcmcControl)
% Write group-constant information to fid
success = true;
fprintf(fid, 'GROUP_Analysis = %d\n', N > 1);
fprintf(fid, 'COVAR_Names = names.txt\n');
fprintf(fid, 'COVAR_Matrix = covar.dat\n');
fprintf(fid, 'SEED_Matrix = seed.dat\n\n');

fprintf(fid, 'POP_Stim = %d\n', size(design, 2));
fprintf(fid, 'Include_temporal_derivative = %d\n\n', ...
    mcmcControl.includeDerivatives);

fprintf(fid, 'MAX_ITER = %d\n', mcmcControl.maxIterations);
fprintf(fid, 'BURN_IN = %d\n', mcmcControl.burnin);
fprintf(fid, 'Expected_Knots = %d\n\n', mcmcControl.expectedKnots);

fprintf(fid, 'NSUBS = %d\n', N);
end


