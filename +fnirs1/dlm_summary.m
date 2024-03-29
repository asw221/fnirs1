
classdef dlm_summary
    % FNIRS1.DLM_SUMMARY stores coefficient summaries for parameters
    % estimated with fnirs1.dlm
    %
    % Example usage:
    %   summ = fnirs1.dlm_summary;
    %   fieldnames(summ);
    %   disp(summ)
    %
    % See also dlm
    % 
    properties
        Nickname;      % Optional, setable, lable for objects
    end
    properties (SetAccess = private)
        Descriptions;  % Description of each parameter
        Estimates;     % Estimate of each parameter
        GroupAnalysis; % Logical group analysis marker
        Intervals;     % Posterior credible intervals for each parameter
        Level;         % Credible interval probability width
        LogFiles;      % List (cellstr) of analysis's MCMC log files
        LogFileDir;    % Analysis directory
        Outcome;       % Outcome type {'hbo', 'hbr', 'hbt'}
        Samples;       % MCMC samples from main analysis
        StdErrors;     % Standard error of each Estimate
    end
    properties (Dependent)
        paramnames;
    end
    properties (Access = private)
        Decimals;            % Max decimal places to display
        DescripPrintWidth;   % Max number of Descriptions characters to print
        EstimatePrintWidth;  % Max numbers to display per item
        FormatDescrip;       % Format specifier for displayed Descriptions
        FormatEstimate;      % Format specifier for displayed Estimates
        FormatEstimateStr;   % Special formatting for displayed Estimates as strings
        FormatInterval;      % Format specifier for displayed Intervals
        FormatLine;          % Format specifier for whole line to display
        IntervalPrintWidth;  % Max width of interval characters to display
        PrintWidth;          % Max number of characters to print for display
        StarsCuttoffs;       % Cuttoff bounds for displaying asterisks
    end
    properties (Hidden, SetAccess = private)
        Burnin;
        Contrasts;     % Set of added contrasts (as a matrix)
    end
    methods
        function obj = dlm_summary()
            % Default constructor for FNIRS1.DLM_SUMMARY objects
            obj.Nickname = '';
            obj.Descriptions = {};
            obj.Estimates = [];
            obj.Intervals = [];
            obj.Level = 0.95;
            obj.LogFiles = {};
            obj.LogFileDir = {};
            obj.Outcome = 'hbo';
            obj.Samples = '';
            obj.StdErrors = [];
            
            obj.Decimals = 3;
            obj.DescripPrintWidth = 45;
            obj.EstimatePrintWidth = 10;
            obj.IntervalPrintWidth = 15;
            obj.FormatDescrip = sprintf('%%-%is', obj.DescripPrintWidth);
            obj.FormatEstimate = sprintf('%%.%if', obj.Decimals);
            obj.FormatEstimateStr = sprintf('%%%is', obj.EstimatePrintWidth);
            obj.FormatInterval = sprintf('(%%.%if, %%.%if)', ...
                obj.Decimals - 1, obj.Decimals - 1);
            obj.FormatLine = sprintf('%%-%is  %%+%is  %%-%is %%%is %%s', ...
                obj.DescripPrintWidth - 2, obj.EstimatePrintWidth - 1, ...
                obj.EstimatePrintWidth - 2, obj.IntervalPrintWidth - 1);  
                % -1's for the spaces
            obj.PrintWidth = 80;  % not counting asterisks
            obj.StarsCuttoffs = [0.95, 0.99, 0.999];
            
            obj.Burnin = 0;
            obj.Contrasts = NaN(0, 0);
        end
        function obj = add_contrast(obj, C)
            con = fnirs1.contrast(C);
            moreFilesNeeded = false;
            modelParamIndex = ~fnirs1.utils.regexpl(obj(1).Descriptions, ...
                '^Contrast: ');
            if (size(con.Vectors, 1) < size(obj(1).Samples, 2))
                con.Vectors = vertcat(con.Vectors, ...
                    zeros(size(obj(1).Samples, 2) - size(con.Vectors, 1), ...
                    size(con.Vectors, 2)));
            elseif (size(con.Vectors, 1) > size(obj(1).Samples, 2))
                if (obj(1).GroupAnalysis)
                    moreFilesNeeded = true;
                else
                    error('Contrast vectors longer than model parameters');
                end
            end
            if (~moreFilesNeeded)
                % --- Everything we need is in obj(:).Samples -------------
                if (obj(1).GroupAnalysis)
                    con.Names = erase(obj(1).Descriptions(modelParamIndex), ...
                        'Population Effect: ');
                else
                    con.Names = obj(1).Descriptions(modelParamIndex);
                end
                for i = 1:numel(obj)
                    BcHat = obj(i).Samples * con.Vectors;
                    p = (1 - obj(i).Level) / 2 * [1 -1] + [0 1];
                    Q = quantile(BcHat, p);
                    if (size(BcHat, 2) > 1)
                        Q = Q';
                    end
                    obj(i).Descriptions = [cellstr("Contrast: " + ...
                        string(con)); obj(i).Descriptions];
                    obj(i).Estimates = [mean(BcHat)'; obj(i).Estimates];
                    obj(i).Intervals = [Q; obj(i).Intervals];
                    obj(i).StdErrors = [std(BcHat)'; ...
                        obj(i).StdErrors];
                    %
                    % New feature: only add contrasts to saved lists in
                    % cases where contrast dim matches obj(i).Samples
                    obj(i).Contrasts = [ obj(i).Contrasts, con.Vectors ];
                end
            else
                % --- More files are needed -------------------------------
                % Need to load subject-specific _beta.log files. This
                % should only be a possibility for group-level analyses
                if (~obj(1).GroupAnalysis)
                    error('DlmSummary:AddContrast:length2', ...
                        'Contrast vector longer than model parameters');
                end
                if (size(con.Vectors, 1) > sum(modelParamIndex))
                    error('DlmSummary:AddContrast:length3', ...
                        'Contrast vector longer than model parameters');
                end
                
                nV = size(con.Vectors, 1);
                groupLevelIndex = fnirs1.utils.regexpl(...
                    obj(1).Descriptions(modelParamIndex), ...
                    '^Population Effect:');
                
                % Find parameters not associated with obj.Samples
                paramsWithMissingSamples = obj(1).Descriptions(...
                    modelParamIndex);
                paramsWithMissingSamples = paramsWithMissingSamples(1:nV);
                groupLevelIndex = groupLevelIndex(1:nV);
                paramsWithMissingSamples = paramsWithMissingSamples(...
                    ~groupLevelIndex);
                
                % Construct list of files we need to look for
                filesNeeded = cell(size(paramsWithMissingSamples));
                for i = 1:numel(filesNeeded)
                    filesNeeded{i} = deblank(strsplit(...
                        paramsWithMissingSamples{i}, 'Cond_'));
                    filesNeeded{i} = sprintf('%s_beta.log', ...
                        filesNeeded{i}{1});
                end
                filesNeeded = unique(filesNeeded);
                
                % Loop over channels, try to find necessary MCMC output
                for i = 1:numel(obj)
                    conCpy = con;
                    try
                        chDir = fileparts(obj(i).LogFileDir);
                        zipDir = sprintf('%s.zip', chDir);
                        chDirExists = exist(chDir, 'dir');
                        zipDirExists = exist(zipDir, 'file');
                        if (zipDirExists && ~chDirExists)
                            unzip(zipDir, fileparts(chDir));
                        end
                        chFilesNeeded = fullfile(obj(i).LogFileDir, ...
                            filesNeeded);
                        datC = cell(size(chFilesNeeded));
                        for j = 1:numel(chFilesNeeded)
                            datC{j} = load(chFilesNeeded{j});
                        end
                        X = [obj(i).Samples, cat(2, datC{:})];
                        conCpy.Vectors = vertcat(con.Vectors, ...
                            zeros(size(X, 2) - nV, size(con.Vectors, 2)));
                        conCpy.Names = obj(i).Descriptions(modelParamIndex);
                        
                        % Add contrasts to model summary
                        BcHat = X * conCpy.Vectors;
                        p = (1 - obj(i).Level) / 2 * [1 -1] + [0 1];
                        Q = quantile(BcHat, p);
                        if (size(BcHat, 2) > 1)
                            Q = Q';
                        end
                        obj(i).Descriptions = [cellstr("Contrast: " + ...
                            string(conCpy)); obj(i).Descriptions];
                        obj(i).Estimates = [mean(BcHat)'; obj(i).Estimates];
                        obj(i).Intervals = [Q; obj(i).Intervals];
                        obj(i).StdErrors = [std(BcHat)'; ...
                            obj(i).StdErrors];
                        
                        % Remove any unzipped files
                        if (zipDirExists && ~chDirExists) % i.e. beforehand
                            remove_folder_and_contents(chDir);
                        end
                    catch ME
                        warning(ME.identifier, '%s', ME.message);
                    end
                end
            end
        end
        function tbl = credint(obj, levels, varargin)
            % Retrieve Bayesian posterior credible intervals for parameter
            % estimates.
            %
            % Usage:
            % credint(obj, [0.8, 0.9, 0.95]) % 80%, 90%, 95% Intervals
            % credint(obj, [0.8, 0.9], 9)    % Intervals for 9th channel
            % credint(obj, [0.8, 0.9], 9, 3) % Show intervals to 3 decimals
            %
            if (nargin <= 3)
                digits = 3;
            else
                digits = max(fix(varargin{2}), 1);
            end
            if (nargin <= 2)
                channel = 1;
            else
                channel = max(fix(varargin{1}), 1);
            end
            levels = max(min(levels, 1), 0);
            qs = interval2quantile(levels);  % fnirs1 private function
            levelname = strsplit(deblank( ...
                sprintf('%.1f%% ', levels * 100)), ' ')';
            samples = obj(channel).Samples;
            pnames = obj(channel).paramnames(1:size(samples, 2));
            if ~isempty(obj(channel).Contrasts)
                con = fnirs1.contrast(obj(channel).Contrasts, pnames);
                samples = [ samples, samples * con.Vectors ];
                pnames = [ pnames; cellstr(con) ];
            end
            [~, ok] = unique(pnames, 'stable');
            p = size(samples, 2);
            cis = cell(p, numel(levels));
            format = sprintf('(%%.%df, %%.%df)', digits, digits);
            for j = 1:size(qs, 2)
                q = quantile(samples, qs(:, j))';
                for i = 1:size(q, 1)
                    cis{i, j} = sprintf(format, q(i, 1), q(i, 2));
                end
            end
            tbl = array2table(cis(ok, :));
            tbl.Properties.VariableNames = levelname;
            tbl.Properties.RowNames = pnames(ok);
        end
        function disp(obj)
            for j = 1:length(obj)
                P = length(obj(j).Estimates);
                if (~isempty(obj(j).Estimates))
                    header = sprintf(obj(j).FormatLine, 'Parameter', 'Estimate', ...
                        'Std.Err.', sprintf('%.1f%% Cred.Int.', obj(j).Level * 100), ' ');
                    fprintf('\tDLM Summary: %s\n', char(obj(j).Nickname));
                    fprintf('%s\n', header);
                    fprintf('%s\n', repmat('-', 1, obj(j).PrintWidth));
                    for i = 1:P
                        descrip = sprintf(obj(j).FormatDescrip, ...
                            fnirs1.utils.abbreviate(obj(j).Descriptions{i}));
                        descrip = descrip(1:min(length(descrip), obj(j).DescripPrintWidth - 2));
                        est = sprintf(obj(j).FormatEstimate, obj(j).Estimates(i));
                        est = est(1:min(length(est), obj(j).EstimatePrintWidth - 1));
                        est = sprintf(obj(j).FormatEstimateStr, est);
                        se = sprintf(obj(j).FormatEstimate, obj(j).StdErrors(i));
                        se = se(1:min(length(se), obj(j).EstimatePrintWidth - 2));
                        intvl = sprintf(obj(j).FormatInterval, obj(j).Intervals(i, :));
                        nstars = sum(obj(j).Level >= obj(j).StarsCuttoffs) * ...
                            (obj(j).Intervals(i, 1) > 0 || obj(j).Intervals(i, 2) < 0);
                        stars = sprintf('%s', repmat('*', 1, nstars));
                        line = sprintf(obj(j).FormatLine, descrip, est, se, intvl, stars);
                        fprintf('%s\n', line);
                    end
                    fprintf('%s\n\n', repmat('-', 1, obj(j).PrintWidth));
                else
                    obj(j).displayEmptyObject();
                end  % if (~isempty(obj(j).Estimate))
            end  % for ob = obj
        end
        function D = get_data(obj, what)
            % Recover stored data from fitted models. Returns a (nested)
            % struct. Parameter 'what' should be a character string with
            % the name of what data to try to recover. Valid options
            % include:
            %   { 'beta',   'delta',  'DLM',    'eta',   'fit',    'HRF',
            %     'knots',  'nknots', 'rawres', 'stdev', 'stdres', 'veta',
            %     'wdelta', 'Xbeta',  'X',      'Y' }
            %
            % This method operates on vectors of fnirs1.dlm_summary
            % objects, but results *may* require a relatively large amount
            % of memory to store
            %
            options = {'beta'; 'delta'; 'DLM'; 'eta'; 'fit'; 'HRF'; ...
                'knots'; 'nknots'; 'rawres'; 'stdev'; 'stdres'; 'veta'; ...
                'wdelta'; 'Xbeta'; 'X'; 'Y'};
            
            if ~any(strcmp(what, options))
                optstr = strcat(options, {', '});
                optstr = strcat(optstr{:});
                error('dlm_summary:GetData:InvalidWhat',...
                    'dlm_summary/get_data: what must be one of:\n[%s]', ...
                    optstr);
            end
            
            N = numel(obj(1).LogFiles);
            D = repmat(struct('name', '', 'data', []), size(obj));
            for i = 1:numel(obj)
                chDir = fileparts(obj(i).LogFileDir);  % channel directory
                zipped = false;
                if ~exist(obj(i).LogFileDir, 'dir')
                    try
                        unzip([chDir '.zip'], fileparts(chDir));
                        zipped = true;
                    catch ME
                        error(ME.message, '%s', ME.identifier);
                    end
                end
                
                % Try to extract relevant data into nested structs
                subD = repmat(struct('what', '', 'data', []), N, 1);
                try
                    for j = 1:N
                        currentFile = erase(obj(i).LogFiles{j}, ...
                            '_beta.log');
                        currentFile = strcat(currentFile, '_', what, '.log');
                        subD(j) = struct('what', ...
                            fnirs1.utils.basename(erase(currentFile, '.log')), ...
                            'data', load(currentFile, '-ascii'));
                    end
                    D(i) = struct('name', obj(i).Nickname, ...
                        'data', subD);
                catch ME
                    error(ME.message, '%s', ME.identifier);
                end
                
                % Clean up unzipped folder if necessary
                if (zipped)
                    remove_folder_and_contents(chDir);
                end
            end
        end
        function value = get.paramnames(obj)
            % Get model parameter names. Returns cellstr
            value = obj(1).Descriptions( ...
                ~fnirs1.utils.regexpl(obj(1).Descriptions, '^Contrast: '));
            value = erase(value, 'Population Effect: ');
        end
        function G = geweke(obj)
            G = table();
            if ~isempty(obj)
                M = size(obj(1).Samples, 2);
                N = numel(obj);
                paramNms = obj.parameter_names();
                paramNms = paramNms(1:M);
                t = NaN(M, N);  % Geweke statistic
                ESS = NaN(M, N);
                Identifier = cell(M, N);
                Parameter = repmat(paramNms, 1, N);
                for i = 1:N
                    t(:, i) = fnirs1.mcmc.geweke(obj(i).Samples);
                    ESS(:, i) = fnirs1.mcmc.ess(obj(i).Samples);
                    Identifier(:, i) = repmat({obj(i).Nickname}, M, 1);
                    % Parameter(:, i) = obj(i).Descriptions(1:M);
                end
                Identifier = reshape(Identifier, N * M, 1);
                Parameter = reshape(Parameter, N * M, 1);
                ESS = reshape(ESS, N * M, 1);
                t = reshape(t, N * M, 1);
                G = table(Identifier, Parameter, ESS, t);
            end
        end
        function H = gewekeplot(obj)
            % Faceted Geweke sequence plots for primary model parameters.
            % Can be used to decide if burnin length is appropriate
            if isempty(obj)
                error('Object contains no data to plot');
            else
                if (numel(obj) > 1)
                    warning('dlm_summary:PlotMultiChannel', ...
                        '%s %s (''%s'')', ...
                        'Object contains data from multiple channels.', ...
                        'Plotting only the first', obj(1).Nickname);
                end
            end
            M = max(floor(sqrt(size(obj(1).Samples, 2))), 1);
            N = ceil(size(obj(1).Samples, 2) / M);
            H = fnirs1.mcmc.gewekeplot(obj(1).Samples);
            nms = obj(1).parameter_names;
            nms = strrep(nms(1:size(obj(1).Samples, 2)), '_', ' ');
            for i = 1:size(obj(1).Samples, 2)
                subplot(M, N, i);
                title(nms{i});
            end
        end
        function obj = head(obj, varargin)
            % Extract only the first N items from an fnris1.dlm_summary
            % object. By default, N = 6
            %
            N = int16(6);
            if (nargin > 1 && isnumeric(varargin{1}))
                N = int16(varargin{1});
                if (N <= 0)
                    error('head: N must be a strictly positive integer');
                end
            end
            for i = 1:numel(obj)
                obj(i).Descriptions = obj(i).Descriptions(1:N);
                obj(i).Estimates = obj(i).Estimates(1:N);
                obj(i).Intervals = obj(i).Intervals(1:N, :);
                obj(i).StdErrors = obj(i).StdErrors(1:N);
            end
        end
        function F = hrf(obj)
            % Return the HRF function (with derivatives) from the fitted
            % model
            F = obj(1).get_data('HRF');
        end
        function B = isempty(obj)
            % Returns logical true of object does not contain any data
            B = isempty(horzcat(obj(:).Estimates));
        end
        function N = parameter_names(obj)
            % Return cellstr of model parameter names
            N = obj.paramnames;
        end
        function obj = parse_log_directory(obj, log_dir)
            % Assign object properties from log files
            
            obj.LogFileDir = fnirs1.utils.explicit_path(log_dir);
            obj.Nickname = fnirs1.utils.basename(fileparts(obj.LogFileDir));
            
            lfd = dir(fullfile(obj.LogFileDir, '*_beta.log'));
            if (isempty(lfd))
                error('dlm_summary:ParseLogs:EmptyLogDir', ...
                    'Empty log file directory (%s)', obj.LogFileDir);
            end
            obj.LogFiles = {lfd(:).name}';
            obj.LogFiles = obj.LogFiles(~fnirs1.utils.regexpl(...
                obj.LogFiles, '^sub_'));
            obj.LogFiles = obj.LogFiles(~fnirs1.utils.regexpl(...
                obj.LogFiles, 'pop_beta.log'));
            obj.LogFiles = fullfile(obj.LogFileDir, obj.LogFiles);
            
            obj.GroupAnalysis = numel(obj.LogFiles) >= 2;
            % Hacky, but should only be 1 *_beta.log file per subj
            obj.Outcome = 'hbo';
            if any(fnirs1.utils.regexpl(obj.LogFiles, '_hbr_'))
                obj.Outcome = 'hbr';
            elseif any(fnirs1.utils.regexpl(obj.LogFiles, '_hbr_'))
                obj.Outcome = 'hbr';
            end
            
            obj.Level = 0.95;
            qs = [(1 - obj.Level) / 2, obj.Level + (1 - obj.Level) / 2];
            
            % Find burnin
            setupfile = fullfile(obj.LogFileDir, '..', 'setup.dat');
            if exist(setupfile, 'file')
                setupinfo = fnirs1.utils.read_setup(setupfile);
                if isfield(setupinfo, 'BURN_IN')
                    obj.Burnin = setupinfo.BURN_IN;
                end
            end
            
            if (obj.GroupAnalysis)
                obj.Samples = load(...
                    fullfile(obj.LogFileDir, 'pop_beta.log'), '-ascii');
                if (obj.Burnin < size(obj.Samples, 1))
                    obj.Samples = obj.Samples((obj.Burnin+1):end, :);
                end
                
                obj.Descriptions = fnirs1.utils.read_whole_file(...
                    fullfile(fileparts(obj.LogFileDir), 'names.txt'));
                obj.Descriptions = cellstr("Population Effect: " + ...
                    string(obj.Descriptions));
                nPopParams = numel(obj.Descriptions);
                
                hasTempDeriv = any(fnirs1.utils.regexpl(...
                    obj.Descriptions, 'TempDeriv'));
                hasDispDeriv = any(fnirs1.utils.regexpl(...
                    obj.Descriptions, 'DispDeriv'));
                
                obj.Estimates = mean(obj.Samples)';
                obj.Intervals = quantile(obj.Samples, qs)';
                obj.StdErrors = std(obj.Samples)';
                
                % Fill out subject level data
                subjBeta = load(obj.LogFiles{1}, '-ascii');
                P = size(subjBeta, 2);
                N = numel(obj.LogFiles);
                subjCoefNames = repmat(" Cond_", P, 1);
                if (hasDispDeriv)
                    subjCoefNames = subjCoefNames + ...
                        repmat(string(1:fix(P/3))', 3, 1) + " HRF";
                    subjCoefNames((fix(P/3) + 1):(2 * fix(P/3))) = ...
                        subjCoefNames((fix(P/3) + 1):(2 * fix(P/3))) + ...
                        " TempDeriv";
                    subjCoefNames((2 * fix(P/3) + 1):end) = ...
                        subjCoefNames((2 * fix(P/3) + 1):end) + ...
                        " DispDeriv";
                elseif (hasTempDeriv)
                    subjCoefNames = subjCoefNames + ...
                        repmat(string(1:fix(P/2))', 2, 1) + " HRF";
                    subjCoefNames((fix(P/2) + 1):end) = ...
                        subjCoefNames((fix(P/2) + 1):end) + " TempDeriv";
                else
                    subjCoefNames = subjCoefNames + ...
                        string(1:P)' + " HRF";
                end
                subjBlockStart = nPopParams + 1;
                obj.Descriptions{N * P + nPopParams} = '';
                obj.Estimates(N * P + nPopParams) = 0;
                obj.Intervals(N * P + nPopParams, :) = 0;
                obj.StdErrors(N * P + nPopParams) = 0;
                for i = 1:N
                    if (i > 1)
                        subjBeta = load(obj.LogFiles{i}, '-ascii');
                    end
                    subjBlockEnd = subjBlockStart + P - 1;
                    bn = erase(fnirs1.utils.basename(obj.LogFiles{i}), ...
                        '_beta.log');
                    obj.Descriptions(subjBlockStart:subjBlockEnd) = ...
                        cellstr(string(bn) + subjCoefNames);
                    obj.Estimates(subjBlockStart:subjBlockEnd) = ...
                        mean(subjBeta)';
                    obj.Intervals(subjBlockStart:subjBlockEnd, :) = ...
                        quantile(subjBeta, qs)';
                    obj.StdErrors(subjBlockStart:subjBlockEnd) = ...
                        std(subjBeta)';
                    subjBlockStart = subjBlockStart + P;
                end
            else
                % Parse single subject model
                obj.Samples = load(obj.LogFiles{1}, '-ascii');
                if (obj.Burnin < size(obj.Samples, 1))
                    obj.Samples = obj.Samples((obj.Burnin+1):end, :);
                end
                
                meanParamDescrip = fnirs1.utils.read_whole_file(...
                    fullfile(fileparts(obj.LogFileDir), 'names.txt'));
                
                hasTempDeriv = any(fnirs1.utils.regexpl(...
                    meanParamDescrip, 'TempDeriv'));
                hasDispDeriv = any(fnirs1.utils.regexpl(...
                    meanParamDescrip, 'DispDeriv'));
                
                obj.Estimates = mean(obj.Samples)';
                obj.Intervals = quantile(obj.Samples, qs)';
                obj.StdErrors = std(obj.Samples)';
                
                % Set descriptions
                P = size(obj.Samples, 2);
                subjCoefNames = repmat(" Cond_", P, 1);
                if (hasDispDeriv)
                    subjCoefNames = subjCoefNames + ...
                        repmat(string(1:fix(P/3))', 3, 1) + " HRF";
                    subjCoefNames((fix(P/3) + 1):(2 * fix(P/3))) = ...
                        subjCoefNames((fix(P/3) + 1):(2 * fix(P/3))) + ...
                        " TempDeriv";
                    subjCoefNames((2 * fix(P/3) + 1):end) = ...
                        subjCoefNames((2 * fix(P/3) + 1):end) + ...
                        " DispDeriv";
                elseif (hasTempDeriv)
                    subjCoefNames = subjCoefNames + ...
                        repmat(string(1:fix(P/2))', 2, 1) + " HRF";
                    subjCoefNames((fix(P/2) + 1):end) = ...
                        subjCoefNames((fix(P/2) + 1):end) + " TempDeriv";
                else
                    subjCoefNames = subjCoefNames + ...
                        string(1:P)' + " HRF";
                end
                bn = erase(fnirs1.utils.basename(obj.LogFiles{1}), ...
                    '_beta.log');
                obj.Descriptions = cellstr(string(bn) + subjCoefNames);
            end
            
            % New add: storage for contrast vectors
            obj.Contrasts = zeros(size(obj.Samples, 2), 0);
        end
        function H = plot(obj)
            % Faceted plot of participants' data and the fitted lines
            if isempty(obj)
                error('Object contains no data to plot');
            else
                if (numel(obj) > 1)
                    warning('dlm_summary:PlotMultiChannel', ...
                        '%s %s (''%s'')', ...
                        'Object contains data from multiple channels.', ...
                        'Plotting only the first', obj(1).Nickname);
                end
            end
            Y = obj(1).get_data('Y');
            fitted = obj(1).get_data('fit');
            M = max(1, floor(sqrt(numel(Y.data))));
            N = ceil(numel(Y.data) / M);
            H = figure;
            % Find range of Y
            yrng = [min(Y.data(1).data), max(Y.data(1).data)];
            if (numel(Y.data) > 1)
                for i = 1:numel(Y.data)
                    yrng = [ min([yrng, Y.data(i).data]), ...
                        max([yrng, Y.data(i).data]) ];
                end
            end
            for i = 1:numel(Y.data)
                subplot(M, N, i);
                tsp = plot(1:numel(Y.data(i).data), Y.data(i).data, ...
                    'LineWidth', 1.6);
                tsp.Color = [0.1, 0.1, 0.1, 0.2];
                hold on
                xlim([1, numel(Y.data(i).data)]);
                ylim(yrng);
                fp = plot(1:numel(fitted.data(i).data), ...
                    fitted.data(i).data, '-m', 'LineWidth', 1.1);
                fp.Color = [0.8, 0, 0.8, 0.6];
                if ( i >= ((numel(Y.data) - N + 1)) )
                    xlabel('Timepoint');
                else
                    set(gca, 'XTickLabel', []);
                end
                if ( rem(i, N) == 1 )
                    ylabel(upper(obj(1).Outcome));
                else
                    set(gca, 'YTickLabel', []);
                end
                title(strrep(fitted.data(i).what, '_', ' '));
            end
        end
        function H = qqplot(obj)
            % Faceted plot of participants' standardized residual against
            % standard Normal quantiles
            if isempty(obj)
                error('Object contains no data to plot');
            else
                if (numel(obj) > 1)
                    warning('dlm_summary:PlotMultiChannel', ...
                        '%s %s (''%s'')', ...
                        'Object contains data from multiple channels.', ...
                        'Plotting only the first', obj(1).Nickname);
                end
            end
            Y = obj(1).get_data('stdres');
            M = max(1, floor(sqrt(numel(Y.data))));
            N = ceil(numel(Y.data) / M);
            H = figure;
            % Find range of Y
            yrng = [min(Y.data(1).data), max(Y.data(1).data)];
            if (numel(Y.data) > 1)
                for i = 1:numel(Y.data)
                    yrng = [ min([yrng, Y.data(i).data]), ...
                        max([yrng, Y.data(i).data]) ];
                end
            end
            for i = 1:numel(Y.data)
                lenY = numel(Y.data(i).data);
                p = (1:lenY) / (lenY + 1);
                z = norminv(p);
                subplot(M, N, i);
                plot(z, sort(Y.data(i).data), 'k.');
                hold on
                plot([z(1), z(end)], [z(1), z(end)], '-r');
                xlim(1.02 * [z(1), z(end)]);
                ylim(1.02 * yrng);
                xlabel('Normal Quantiles');
                ylabel('Residual Quantiles');
                title(strrep(Y.data(i).what, '_', ' '));
            end
        end
        function obj = read_from_file(obj, file)
            % Read data from a Parameter_Estimates.log file and extract
            % parameter estimates and inferential summaries
            
            warning('dlm_summary read_from_file method is deprecated');
            
            obj.LogFileDir = fileparts(file);
            if (isempty(obj.LogFileDir) || ...
                    strcmp(obj.LogFileDir, './') || ...
                    strcmp(obj.LogFileDir, file))
                obj.LogFileDir = pwd;
            end
            lfd = dir(fullfile(obj.LogFileDir, '*_beta.log'));
            obj.LogFiles = fullfile(obj.LogFileDir, {lfd(:).name}');
            obj.Nickname = fnirs1.utils.basename(fileparts(obj.LogFileDir));
            
            % Identify outcome type
            outcomeTypes = {};
            if (any(fnirs1.utils.regexpl(obj.LogFiles, '_hbr_')))
                outcomeTypes = [outcomeTypes, {'hbr'}];
                obj.Outcome = 'hbr';
            elseif (any(fnirs1.utils.regexpl(obj.LogFiles, '_hbt_')))
                outcomeTypes = [outcomeTypes, {'hbt'}];
                obj.Outcome = 'hbt';
            elseif (any(fnirs1.utils.regexpl(obj.LogFiles, '_hbo_')))
                outcomeTypes = [outcomeTypes, {'hbo'}];
                obj.Outcome = 'hbo';
            end
            if (isempty(outcomeTypes))
                warning('Outcome type not identified in log files');
            elseif (numel(outcomeTypes) > 1)
                warning('Multiple outcome types appear in log file names: %s', ...
                    strjoin(outcomeTypes));
            end
            
            
            % Read file into cell array and clean up lines
            lines = fnirs1.utils.read_whole_file(file);
            lines = lines(~strcmp(lines, ''));
            if (length(lines) > 2)
                % last two lines of Parameter_Estimates.log files are notes
                % explaining quantile/interval methods
                lines = lines(1:(end - 2));
            end
            
            % Clean lines, increment stimulus/covariate numbers (originally
            % start from zero)
            for i = 1:length(lines)
                lines{i} = erase(lines{i}, ': Parameter Summary for regression of');
                lines{i} = strrep(lines{i}, 'Temporal Derivative', 'TempDeriv');
                lines{i} = strrep(lines{i}, 'Pop ', 'Population ');
                lines{i} = strrep(lines{i}, ...
                    'Population Level Parameters:', 'Population Effect');
                
                condNo = sscanf(lines{i}, 'Cond_%i');
                if ~isempty(condNo)
                    lines{i} = strrep(lines{i}, sprintf('Cond_%i', condNo), ...
                        sprintf('Cond_%i', condNo));
                end
                
                condNo = sscanf(lines{i}, 'Cond = %i');
                if ~isempty(condNo)
                    lines{i} = strrep(lines{i}, sprintf('Cond = %i', condNo), ...
                        sprintf('Cond_%i', condNo + 1));
                end
                
                covarNo = sscanf(lines{i}, '%*s Covariate %i');
                if ~isempty(covarNo)
                    lines{i} = strrep(lines{i}, sprintf('Covariate %i', covarNo), ...
                        sprintf('Covar %i', covarNo + 1));
                end
            end
            
            % Identify special types of lines/blocks related to parameter
            % interpretation - group vs subject-specific analyses
            mean_lines = contains(lines, 'mean = ');
            P = sum(mean_lines);
            if (P == 0)
                error('%s does not appear to be a valid parameter-log file', file);
            end
            
            % File is structured:
            %   Title
            %      SubTitle
            %        mean = X\tsd = Y
            %        V% Cred.Int. = (L, U)*
            interval_lines = [false; mean_lines(1:(end - 1))];
            subtitle_lines = [mean_lines(2:end); false];
            
            % Set object parameter values
            obj.Descriptions = cell(P, 1);
            obj.Estimates = nan(P, 1);
            obj.StdErrors = nan(P, 1);
            obj.Intervals = nan(P, 2);
            
            count = 1;
            title = '';
            blockTitle = '';
            groupAnalysisBlock = true(length(lines), 1);
            beginningOfSubjectBlockLine = 'Parameter Summary for Subjects';
            beginningOfSubjectBlockIndex = find(contains(lines, ...
                beginningOfSubjectBlockLine));
            if (~isempty(beginningOfSubjectBlockIndex))
                groupAnalysisBlock(beginningOfSubjectBlockIndex(1):end) = false;
            end
            obj.GroupAnalysis = any(groupAnalysisBlock);
            
            subtitle_lines(groupAnalysisBlock) = false;
            title_lines = ~(mean_lines | interval_lines | subtitle_lines);
            
            % Loop over lines and extract relevant information.
            % (!) This block may need to be edited if Parameter_Estimates.log
            % formatting ever chanes
            for i = 1:length(lines)
                % Identify line types
                if (title_lines(i))
                    % Identify title/subtitle information - formatting slightly
                    % different for group analysis results vs subject-specific
                    % analyses
                    if (strcmpi(lines{i}, beginningOfSubjectBlockLine))
                        title_lines(i) = false;
                    end
                    if (i == 1)
                        blockTitle = lines{i};
                    else
                        if (length(lines) > i)
                            if (~groupAnalysisBlock(i) && ...
                                    interval_lines(i - 1) && ...
                                    title_lines(i + 1))
                                blockTitle = lines{i};
                            elseif (groupAnalysisBlock(i) && ...
                                    interval_lines(i - 1) || ...
                                    strcmpi(lines{i - 1}, ...
                                    beginningOfSubjectBlockLine))
                                blockTitle = lines{i};
                            else
                                title = lines{i};
                            end
                        end
                    end
                    if ~isempty(blockTitle)
                        blockTitle = strrep(blockTitle, ':', '');
                    end
                elseif (subtitle_lines(i))
                    % Format parameter descriptions given block-titles,
                    % titles, and subtitles
                    obj.Descriptions{count} = strrep(...
                        sprintf('%s %s: %s', blockTitle, title, lines{i}), ...
                        ' :', ':');
                    if (contains(obj.Descriptions{count}, 'Population, '))
                        obj.Descriptions{count} = sprintf('Population %s', ...
                            erase(obj.Descriptions{count}, 'Population, '));
                    end
                elseif (mean_lines(i))
                    theta = format_mean_sd_line(lines{i});
                    obj.Estimates(count) = theta.Estimate;
                    obj.StdErrors(count) = theta.SE;
                    if (~isempty(theta.Name))
                        obj.Descriptions{count} = strrep(...
                            sprintf('%s %s: %s', blockTitle, title, theta.Name), ...
                            ' :', ':');
                        if (contains(obj.Descriptions{count}, 'Population, '))
                            obj.Descriptions{count} = sprintf('Population %s', ...
                                erase(obj.Descriptions{count}, 'Population, '));
                        end
                    end
                elseif (interval_lines(i))
                    ntrvl = format_interval_line(lines{i});
                    obj.Intervals(count, :) = ntrvl.Interval';
                    obj.Level = ntrvl.Level;
                    count = count + 1;
                end
            end  % end loop over lines
            
            % If TempDeriv's included, awkwardness that MCMC samples
            % columns are in order: [HRF_1-n, TempDeriv_1-n], but
            % Parameter_Estimates.log order is interleaved: [HRF_1,
            % TempDeriv_1, ..., HRF_n, TempDeriv_n]. Reorder Estimates,
            % SE's and intervals:
            derivInd = fnirs1.utils.regexpl(...
                obj.Descriptions, 'Deriv') & ...
                ~fnirs1.utils.regexpl(obj.Descriptions, 'Population');
            if (any(derivInd))
                obj.Descriptions = [obj.Descriptions(~derivInd); ...
                    obj.Descriptions(derivInd)];
                obj.Estimates = [obj.Estimates(~derivInd); ...
                    obj.Estimates(derivInd)];
                obj.StdErrors = [obj.StdErrors(~derivInd); ...
                    obj.StdErrors(derivInd)];
                obj.Intervals = [obj.Intervals(~derivInd, :); ...
                    obj.Intervals(derivInd, :)];
            end
            
            % Read in appropriate samples
            if (obj.GroupAnalysis)
                try
                    obj.Samples = load(...
                        fullfile(obj.LogFileDir, 'pop_beta.log'), ...
                        '-ascii');
                catch ME
                    warning(ME.identifier, '%s', ME.message);
                end
            else
                % Attempt to find correct file for single subject
                obj.Samples = 0;
                for i = 1:length(lfd)
                    if ~(fnirs1.utils.regexpl(lfd(i).name, '^sub_') || ...
                            fnirs1.utils.regexpl(lfd(i).name, '^pop_beta'))
                        if all(obj.Samples == 0)
                            obj.Samples = load(obj.LogFiles{i}, '-ascii');
                        end
                    end
                end
                if all(obj.Samples == 0)
                    obj.Samples = [];
                    warning('Could not locate participant''s MCMC ouptut file');
                end
            end
        end
        function H = residplot(obj)
            % Faceted plot of participants' residual time series
            if isempty(obj)
                error('Object contains no data to plot');
            else
                if (numel(obj) > 1)
                    warning('dlm_summary:PlotMultiChannel', ...
                        '%s %s (''%s'')', ...
                        'Object contains data from multiple channels.', ...
                        'Plotting only the first', obj(1).Nickname);
                end
            end
            Y = obj(1).get_data('rawres');
            M = max(1, floor(sqrt(numel(Y.data))));
            N = ceil(numel(Y.data) / M);
            H = figure;
            % Find range of Y
            yrng = [min(Y.data(1).data), max(Y.data(1).data)];
            if (numel(Y.data) > 1)
                for i = 1:numel(Y.data)
                    yrng = [ min([yrng, Y.data(i).data]), ...
                        max([yrng, Y.data(i).data]) ];
                end
            end
            for i = 1:numel(Y.data)
                subplot(M, N, i);
                plot(1:numel(Y.data(i).data), Y.data(i).data, 'k.');
                hold on
                plot([1, numel(Y.data(i).data)], [0, 0], '-r');
                xlim([1, numel(Y.data(i).data)]);
                ylim(yrng);
                xlabel('Timepoint');
                ylabel([upper(obj(1).Outcome) ' Residual']);
                title(strrep(Y.data(i).what, '_', ' '));
            end
        end
        function obj = set.Nickname(obj, nname)
            try
                nname = char(nname);
                obj.Nickname = nname;
            catch ME
                warning(ME.identifier, '%s', ME.message);
                warning('DlmSummary:SetNickname', ...
                    'Nickname property must be coercible to char');
            end
        end
        function tbl = table(obj, varargin)
            % Conversion from FNIRS1.DLM_SUMMARY to table
            %
            % Optional additional argument can be used to select whether to
            % return table of the estimates or standard errors. For
            % example, if D is an fnirs1.dlm_summary object,
            % >> table(D)              % returns table of Estimates
            % >> table(D, 'StdError')  % returns table of standard errors
            % >> table(D, 'ZStat')     % returns table of z statistics
            %
            prop = 'Estimates';
            if (nargin > 1)
                if ~(ischar(varargin{1}) || isstring(varargin{1}))
                    error('dlm_summary:TableConversion:BadProperty', ...
                        'property should be string-like');
                end
                if (strcmpi(varargin{1}, 'StdError') || ...
                        strcmpi(varargin{1}, 'Std.Err.') || ...
                        strcmpi(varargin{1}, 'StdErrors') || ...
                        strcmpi(varargin{1}, 'SE'))
                    prop = 'StdErrors';
                elseif (strcmpi(varargin{1}, 'ZStat') || ...
                        strcmpi(varargin{1}, 'Z') || ...
                        strcmpi(varargin{1}, 'TStat') || ...
                        strcmpi(varargin{1}, 'T') || ...
                        strcmpi(varargin{1}, 'Intensity'))
                    prop = 'ZStat';
                elseif ~(strcmpi(varargin{1}, 'Estimate') || ...
                        strcmpi(varargin{1}, 'Estimates'))
                    error('dlm_summary:TableConversion:UnrecProp', ...
                        'unrecognized property option');
                end
            end
            if ~isempty(obj)
                Channel = vertcat(obj(:).Nickname);
                if (isempty(Channel))
                    Channel = string(1:numel(obj));
                end
                Channel = cellstr(Channel);
                if (strcmp(prop, 'ZStat'))
                    tblData = horzcat(obj(:).Estimates)' ./ ...
                        horzcat(obj(:).StdErrors)';
                else
                    tblData = horzcat(obj(:).(prop))';
                end
                tbl = array2table(tblData, ...
                    'VariableNames', ...
                    matlab.lang.makeValidName(obj(1).Descriptions), ...
                    'RowNames', Channel);
            else
                tbl = table();
            end
        end
        function H = traceplot(obj)
            % Create MCMC traceplots of primary model parameters
            if isempty(obj)
                error('Object contains no data to plot');
            else
                if (numel(obj) > 1)
                    warning('dlm_summary:PlotMultiChannel', ...
                        '%s %s (''%s'')', ...
                        'Object contains data from multiple channels.', ...
                        'Plotting only the first', obj(1).Nickname);
                end
            end
            M = max(floor(sqrt(size(obj(1).Samples, 2))), 1);
            N = ceil(size(obj(1).Samples, 2) / M);
            H = figure;
            nms = strrep(obj(1).parameter_names(), '_', ' ');
            T = size(obj(1).Samples, 1);
            for i = 1:size(obj(1).Samples, 2)
                subplot(M, N, i);
                plot(1:T, obj(1).Samples(:, i), '-k');
                xlabel('Iteration');
                title(nms{i});
            end
        end
    end
    methods (Static)
        function displayEmptyObject()
            fprintf('\tempty fnirs1.dlm_summary object\n');
        end
    end
end
