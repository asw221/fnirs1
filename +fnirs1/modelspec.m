
classdef modelspec
    % Primary interface to create and manage model prerequisite
    % information. Users can also create fnirs1.modelspec objects through
    % the fnirs1.specify_model interface, which provides slightly more
    % flexibility.
    %
    % Example usage:
    %   spec = fnirs1.modelspec(datafiles);  % datafiles lists *.mat files
    %   spec = spec.write_channel_data() ;   % Reassignment is important
    %   fit = fnirs1.dlm(spec);
    %
    % See also fnirs1.specify_model
    %
    properties
        channels;   % Indices of channel subsets
        credints;   % Credible interval widths
        downsamp;   % Rate data should be down-sampled
        mcmc;       % MCMC control
        outcome;    % Outcome type: one of {'hbo', 'hbt', 'hbr'}
    end
    properties (SetAccess = private)
        datafiles;  % Cell array of outcome file names
        gppar;      % Parameters for group analysis (empty for non-group)
        nchan;      % Number of channels collected
        nstim;      % Number of task conditions/stimulus types
        setupfiles; % Paths to 'setup.dat' files
    end
    properties (Access = private)
        spath;  % Model(s) setup file path
    end
    properties (Dependent)
        group;  % Group analysis boolean
        n;      % Number of participants
        nrep;   % Number of repeats per subject
    end
    methods
        function obj = modelspec(varargin)
            % Default values
            obj.datafiles = {''};
            obj.downsamp = 1;
            obj.channels = [];
            obj.credints = [0.8, 0.9, 0.95, 0.99];
            obj.gppar = struct();
            obj.mcmc = fnirs1.mcmc_control;
            obj.nstim = [];
            obj.outcome = 'hbo';
            if (nargin >= 1)
                obj = obj.parse_arguments(varargin);
                [ns, nch] = obj.get_nstim_nchan();
                obj.nstim = ns;
                obj.nchan = nch;
                obj.spath = obj.create_output_dir();
                if (numel(obj.downsamp) ~= obj.n)
                    obj.downsamp = repmat(obj.downsamp, ...
                        ceil(obj.n / numel(obj.downsamp)), 1);
                    obj.downsamp = obj.downsamp(1:obj.n);
                end
                if isempty(obj.channels)
                    obj.channels = 1:obj.nchan;
                end
                if (obj.group)
                    if (size(obj.gppar.data, 1) ~= numel(obj.datafiles))
                        error('modelspec:mismatcheddata', ...
                            [num2str(numel(obj.datafiles)), ...
                            ' data files provided, but ''GroupData'' ', ...
                            'table has ', ...
                            num2str(size(obj.gppar.data), 1), ...
                            ' rows']);
                    end
                    obj.gppar.data = fnirs1.expand_table_conditions( ...
                        obj.gppar.data, obj.nstim, 'Cond');
                end
            end
        end
        function obj = parse_arguments(obj, args)
            % Parse input name/value pair arguments
            options = struct( ...
                'CredibleIntervals', [], ...
                'DownSampleRate', [], ...
                'GroupCovariateNames', '', ...
                'GroupData', [], ...
                'GroupFormula', '', ...
                'McmcControl', table(), ...
                'OutcomeType', '', ...
                'SpecificChannels', [] ...
                );
            optnames = fieldnames(options);
            obj.datafiles = obj.validate_datafiles(args{1});
            nargs = numel(args) - 1;
            if (mod(nargs, 2) ~= 0)
                error('modelspec:parse_arguments:unpairedproperty', ...
                'Pair/value mismatch');
            end
            for pair = reshape(args(2:end), 2, [])
                nm = pair{1};
                if (any(strcmp(nm, optnames)))
                    options.(nm) = pair{2};
                else
                    error('modelspec:parse_arguments:unknownargument', ...
                        '%s is not a recognized parameter', nm);
                end
            end
            % set any requested parameters
            if ~isempty(options.CredibleIntervals)
                obj.credints = options.CredibleIntervals;
            end
            if ~isempty(options.DownSampleRate)
                obj.downsamp = options.DownSampleRate;
            end
            if ~isempty(options.McmcControl)
                obj.mcmc = options.McmcControl;
            end
            if ~isempty(options.OutcomeType)
                obj.outcome = options.OutcomeType;
            end
            if ~isempty(options.SpecificChannels)
                obj.channels = options.SpecificChannels;
            end
            if (obj.group)
                obj.gppar = obj.empty_group_parameters();
                % For now, treat empty data/formula as an error
                if (isempty(options.GroupData) || ...
                        isempty(options.GroupFormula))
                    error(['Cannot have missing ''GroupData'' or ', ...
                        '''GroupFormula'' for group analyses.']);
                end
                %
                if ~isempty(options.GroupCovariateNames)
                    warning('Setting ''GroupCovariateNames'' is deprecated');
                    % obj.gppar.covnames = options.GroupCovariateNames;
                end
                if ~isempty(options.GroupData)
                    if isempty(options.GroupFormula)
                        error('modelspec:parse_arguments:missingformula', ...
                            ['''GroupData'' parameter must be ', ...
                            'accompanied by a ''GroupFormula'' parameter']);
                    end
                    obj.gppar.data = options.GroupData;
                end
                if ~isempty(options.GroupFormula)
                    if isempty(options.GroupData)
                        error('modelspec:parse_arguments:missingdata', ...
                            ['''GroupFormula'' parameter must be ', ...
                            'accompanied by a ''GroupData'' parameter']);
                    end
                    obj.gppar.formula = options.GroupFormula;
                end
            end
        end
        function [nstim, nchan] = get_nstim_nchan(obj)
            try
                data = load(obj.datafiles{1});
                nstim = size(data.s, 2);
                nchan = size(data.(obj.outcome), 2);
            catch me
                error('modelspec:get_nstim_nchan', ...
                    ['Could not deduce condition/channel numbers ', ...
                    'from input\n\t-> (%s)'], me.message);
            end
        end
        function value = get.group(obj)
            value = size(obj.datafiles, 1) > 1;
        end
        function value = get.n(obj)
            value = size(obj.datafiles, 1);
        end
        function value = get.nrep(obj)
            value = size(obj.datafiles, 2);
        end
        function value = group_design(obj, varargin)
            % Return a struct with fields 'x'---the group design matrix,
            % and 'names'---the names of the variables in the columns of
            % 'x'.
            % Optional argument to indicate if covariates should be put in
            % {Cond_*, TempDeriv_*, DispDeriv_*, Interaction} order
            % (default = true)
            reorder = true;
            if (nargin > 1)
                if islogical(varargin{1})
                    reorder = varargin{1};
                else
                    warning('modelspec:group_design:ignoredarg', ...
                        'Ignored argument (unrecognized type)');
                end
            end
            %
            d = fnirs1.mixed_effects_design(obj.gppar.data, ...
                obj.gppar.formula);
            if (fnirs1.utils.regexpl(obj.gppar.formula, 'Cond'))
                d = d.reformat_base_condition('Cond');
            end
            value = struct('x', designMatrix(d, 'Fixed'), ...
                'names', string(d.CoefficientNames));
            %
            % Append HRF derivatives if requested
            if (obj.mcmc.hrfDerivatives > 0)
                value = obj.append_hrf_derivatives(value);
            end
            if (reorder)
                value = obj.reorder_group_design(value);
            end
        end
        function obj = init_setups_and_channeldirs(obj)
            % Initialize output directories with setup.dat's
            obj.setupfiles = cell(numel(obj.channels), 1);
            for c = 1:obj.channels
                chdir = fullfile(obj.spath, ...
                    sprintf('ch%04d', obj.channels(c)));
                obj.setupfiles{c} = fullfile(chdir, 'setup.dat');
                [~, ~, ~] = mkdir(chdir);
                % ^^ strange syntax suppresses warnings if 'chdir' already
                % exists
                obj.write_setup_preamble(obj.setupfiles{c});
                write_seed(chdir);  % fnirs1 private function
            end
            %
            % Write covariate names and group covariate matrix
            if (obj.group)
                design = obj.group_design();
                for c = 1:numel(obj.setupfiles)
                    chdir = fileparts(obj.setupfiles{c});
                    obj.write_group_data(chdir, design);
                end
                % obj.write_group_data(chdir);
            end
        end
        function obj = set.channels(obj, value)
            if all(value >= 0)
                obj.channels = int32(value);
            else
                error('modelspec:set:channels', ...
                    'Cannot set channel subset to indices < 1');
            end
        end
        function obj = set.credints(obj, value)
            obj.credints = max(min(value, 1), 0);
        end
        function obj = set.downsamp(obj, value)
            if all(value > 0)
                obj.downsamp = int32(value);
            else
                error('modelspec:set:downsamp', ...
                    'Cannot set negative downsample rate');
            end
        end
        function obj = set.mcmc(obj, value)
            if isa(value, 'fnirs1.mcmc_control')
                obj.mcmc = value;
            else
                error('modelspec:mcmc', ...
                    'mcmc must be an fnirs1.mcmc_control object');
            end
        end
        function obj = set.outcome(obj, value)
            options = {'hbo', 'hbt', 'hbr'};
            value = cellstr(value);
            value = char(value{1});
            if any(strcmp(options, value))
                obj.outcome = value;
            else
                error('modelspec:outcome:unrecognizedoption', ...
                    'Unrecognized outcome: %s', value);
            end
        end
        function obj = write_channel_data(obj)
            % Make sub-directories for each channel and write the data
            obj = obj.init_setups_and_channeldirs();
            warnings_generated = false;
            for i = 1:obj.n
                tasknames = "task" + string(1:obj.nrep);
                pdatafiles = "";
                pstimfiles = "";
                for r = 1:obj.nrep
                    % Save participant's data file basename
                    [~, bn] = fileparts(obj.datafiles{i, r});
                    % Make sure data can be loaded and that the file has
                    % the correct fields
                    try data = load(obj.datafiles{i, r});
                        if (size(data.s, 2) ~= obj.nstim)
                            error(['Mismatch in number of task conditions: ', ...
                                '''%s'' and ''%s''. fnirs1 does not support ', ...
                                'imbalanced designs'], ...
                                obj.datafiles{1, 1}, obj.datafiles{i, r});
                        end
                        if obj.missing_fields(data, obj.outcome)
                            [~, msng] = obj.missing_fields(data, obj.outcome);
                            error('File ''%s''\n\tis missing fields:  %s\b\b', ...
                                obj.datafiles{i, r}, ...
                                sprintf('%s, ', string(msng)));
                        end
                        if (size(data.(obj.outcome), 1) ~= size(data.s, 1))
                            error('modelspec:write_channel_data:timestamp', ...
                                ['File %s timepoint mismatch:\n\tField ', ...
                                '''%s'' has %d entries while ''s'' has %d'], ...
                                fnirs1.utils.basename(obj.datafiles{i, r}), ...
                                obj.outcome, ...
                                size(data.(obj.outcome), 1), ...
                                size(data.s, 1) );
                        end
                        if (any(obj.channels) > size(data.(obj.outcome), 2))
                            error('modelspec:write_channel_data:nchannels', ...
                                ['File %s has %d channels recorded', ...
                                '\b\t-> (requested channels: %s)'], ...
                                obj.datafiles{i, r}, ...
                                size(data.(obj.outcome), 2), ...
                                sprintf('%d, ', obj.channels));
                        end
                    catch me
                        % remove_folder_and_contents(obj.spath);
                        error(me.identifier, '%s', me.message);
                    end
                    samplerate = obj.fnirs_sampling_rate(data);
                    dsr = obj.downsamp(i);
                    if ~obj.valid_downsample_rate(samplerate, dsr)
                        error('modelspec:write_channel_data:downsamplerate', ...
                            ['Improper down-sampling rate (%.2f) for ', ...
                            'data from file:\n\t%s'], ...
                            dsr, obj.datafiles{i, r});
                    end
                    % Write subject specific channel data
                    y = data.(obj.outcome);
                    pdatafiles(r) = string(bn) + "_" + ...
                        tasknames(r) + "_" + string(obj.outcome) + ".txt";
                    pstimfiles(r) = string(bn) + "_" + ...
                        tasknames(r) + "_stim.txt";
                    for c = 1:numel(obj.channels)
                        chan = obj.channels(c);
                        chdir = fileparts(obj.setupfiles{c});
                        pdf = fullfile(chdir, pdatafiles(r));
                        psf = fullfile(chdir, pstimfiles(r));
                        if (~write_matrix(y(:, chan), pdf) || ...
                                ~write_matrix(data.s, psf))
                            % write_matrix is another fnirs1 private func
                            warning('modelspec:writedatafiles', ...
                                ['Error writing data or stim file ', ...
                                'for: ', bn]);
                            warnings_generated = true;
                        end
                        %
                        % Append participant info to [channel]/setup.dat
                        % files
                        if (r == obj.nrep)
                            obj.write_setup_participant( ...
                                obj.setupfiles{c}, ...
                                pdatafiles, pstimfiles, samplerate, dsr);
                        end
                    end  % for c = 1:numel(obj.channels)
                    fprintf('\tFiles written for %s (task %d)\n', bn, r);
                end  % for r = 1:obj.nrep
            end  % for i = 1:obj.n
            if (warnings_generated)
                % remove_folder_and_contents(obj.spath);
                error(['Warnings generated. Please double check ', ...
                    'participant data. Setup files have NOT been removed.']);
            end
        end
        function write_setup_preamble(obj, fname)
            fid = fopen(fname, 'w');
            if (fid < 0)
                error('modelspec:write_setup_preamble:invalidfid', ...
                    'Setup file not open');
            end
            fprintf(fid, 'GROUP_Analysis = %.f\n', obj.group);
            fprintf(fid, 'COVAR_Names = names.txt\n');
            fprintf(fid, 'COVAR_Matrix = covar.dat\n');
            fprintf(fid, 'SEED_Matrix = seed.dat\n\n');
            
            fprintf(fid, 'POP_Stim = %.f\n', obj.nstim);
            fprintf(fid, 'Include_temporal_derivative = %.f\n', ...
                obj.mcmc.hrfDerivatives);
            fprintf(fid, 'Confidence_Intervals = %s\n\n', ...
                sprintf('%.6f ', obj.credints));
            
            fprintf(fid, 'MAX_ITER = %.f\n', obj.mcmc.maxIterations);
            fprintf(fid, 'BURN_IN = %.f\n', obj.mcmc.burnin);
            fprintf(fid, 'Expected_Knots = %.f\n\n', obj.mcmc.expectedKnots);
            
            fprintf(fid, 'NSUBS = %.f\n', obj.n);
            % fprintf(obj.sid, 'Random = %.f\n', opts.NRandom);
            % fprintf(obj.sid, 'Fixed = %.f\n', opts.NFixed);
            fclose(fid);
        end
    end
    methods (Static)
        function p = create_output_dir()
            % Create (temporary?) directory for output. Directory will have
            % a complex name to avoid potential clashes if "parallelizing"
            % by hand.
            p = fullfile(pwd, sprintf('_fnirs_%s_%s', ...
                datetime('now', 'TimeZone', 'local', 'Format', ...
                'd-MMM-y-HH-mm-ss-SSS'), ...
                fnirs1.utils.basename(tempname)) );
            [~, ~, ~] = mkdir(p);
            % ^^ strange syntax suppresses warnings if 'p' already exists
        end
        function s = empty_group_parameters()
            s = struct('data', table(), 'formula', '');
        end
        function f = parse_group_formula(formula)
            % remove any instance of TempDeriv from the formula
            if (contains(formula, 'TempDeriv'))
                % regular expression is not foolproof (fails on 
                % tertiary+ interactions)
                fprintf('%s -> ', formula);
                formula = regexprep(formula, ...
                    '[ *+-(:]*TempDeriv[):]*', '');
                fprintf('%s\n', formula);
                warning('Don''t put TempDeriv directly in model formula');
            end
            f = formula;
        end
        function r = fnirs_sampling_rate(data)
            if (numel(data.t) < 2)
                error('modelspec:fnirs_sampling_rate', ...
                    'Invalid timestamps');
            end
            r = floor(1 / (data.t(2) - data.t(1)));
        end
        function [tf, names] = missing_fields(data, outcome)
            % Check if any fields are missing from data format. Return true
            % if any fields are missing. As a second return, return the
            % names of the missing fields
            fn = fieldnames(data);
            needed = {outcome; 's'; 't'};
            ok = logical(size(needed));
            for i = 1:numel(needed)
                ok(i) = any(strcmp(needed{i}, fn));
            end
            tf = any(~ok);
            names = needed(~ok);
        end
        function tf = valid_downsample_rate(samplerate, downsamplerate)
            multiple = fix(samplerate / downsamplerate) == ...
                (samplerate / downsamplerate);
            tf = (downsamplerate > 0) && multiple;                
        end
        function f = validate_datafiles(files)
            % Check that input files exist and argument can be converted to
            % cellstr type
            try
                f = cellstr(files);
            catch
                error('modelspec:validate_datafiles:argtype', ...
                    'Input must be cellstr like');
            end
            valid = false(numel(f));
            for i = 1:numel(f)
                valid(i) = logical(exist(f{i}, 'file'));
            end
            if ~all(valid)
                error('modelspec:validate_datafiles:nonexistantfiles', ...
                    'The following files were not found:\n%s\n', ...
                    sprintf('%s\n', string(f(~valid))));
            end
            if (size(f, 1) == 1 && size(f, 2) ~= 1)
                warning('modelspec:validate_datafiles:inputdimension', ...
                    ['Transposing input data files to avoid ', ...
                    'treating analysis as repeated measures on ', ...
                    'a single participant']);
                f = f';
            end
        end
        function write_setup_participant(fname, datafiles, ...
                designfiles, samplingrate, downsamplerate)
            fid = fopen(fname, 'a');
            if (fid < 0)
                error('modelspec:write_setup_participant:invalidfid', ...
                    'Setup file not open');
            end
            fprintf(fid, '\nSUB_Replicates = %d\n', length(datafiles));
            fprintf(fid, 'SUB_Data = %s\n', ...
                deblank(sprintf('%s ', ...
                string(fnirs1.utils.basename(datafiles))) ));
            fprintf(fid, 'SUB_Design = %s\n', ...
                deblank( sprintf('%s ', ...
                string(fnirs1.utils.basename(designfiles))) ));
            fprintf(fid, 'SUB_Freq = %d\n', samplingrate);
            fprintf(fid, 'SubSamp_Freq = %d\n', downsamplerate);
            fclose(fid);
        end
    end
    methods (Hidden)
        function design = append_hrf_derivatives(obj, design)
            % Appends HRF derivatives to group design.
            % Input 'design' should be a struct with fields 'x' and 'names'
            % as output by modelspec::group_design
            if (obj.mcmc.hrfDerivatives > 0)
                nderiv = obj.mcmc.hrfDerivatives * obj.nstim;
                % Append the group design matrix by creating a smaller
                % block diagonal matrix ([X 0; 0 I]) and indexing the rows
                % to insert/repeat the 'I' rows
                xi = blkdiag(design.x, eye(nderiv));
                xrows = ( 1:size(design.x, 1) )';
                irows = ( (size(design.x, 1) + 1):size(xi, 1) )';
                m = fix( size(design.x, 1) / obj.nstim );
                ind = reshape( [ ...
                    reshape(xrows, obj.nstim, m); ...
                    repmat(irows, 1, m) ...
                    ], [], 1);
                design.x = xi(ind, :);
                % append names:
                switch obj.mcmc.hrfDerivatives
                    case 1
                        design.names = [design.names, ...
                            cellstr("TempDeriv_" + string(1:obj.nstim))];
                    case 2
                        design.names = [design.names, ...
                            cellstr("TempDeriv_" + string(1:obj.nstim)), ...
                            cellstr("DispDeriv_" + string(1:obj.nstim))];
                end
            end
        end
    end
    methods (Static, Hidden)
        function design = reorder_group_design(design)
            % Reorder the disign to put the random effects 
            % {Cond_*, TempDeriv_*, DispDeriv_*, Interactions}
            % first
            lcond  = fnirs1.utils.regexpl(design.names, 'Cond_');
            ltempd = fnirs1.utils.regexpl(design.names, 'TempDeriv_');
            ldispd = fnirs1.utils.regexpl(design.names, 'DispDeriv_');
            linter = fnirs1.utils.regexpl(design.names, ':');
            %
            ind = unique( [ ...
                find(lcond  & ~linter), ... % Put ordered Cond_* first
                find(ltempd & ~linter), ... % Then TempDeriv_*
                find(ldispd & ~linter), ... % Then DispDeriv_*
                find(lcond  &  linter), ... % Then interactions with above
                find(ltempd &  linter), ...
                find(ldispd &  linter)  ...
                ], 'stable' );
            % ^^ 'stable' does not alter order
            % Legacy: (Not currently needed)
            % ---
            % reNdx = ind;
            % fixNdx = find(~(lcond | ltempd | ldispd));
            % NRandom = numel(reNdx);
            % NFixed = numel(fixNdx);
            % ---
            ord = [ind, setdiff(1:size(design.x, 2), ind)];
            % Reorder design/names
            design.x = design.x(:, ord);
            design.names = design.names(ord);
        end
        function write_group_data(outdir, design)
            % Write group-level design matrix and covariate names to
            % 'covar.dat' and 'names.txt' files in outdir.
            % Input 'design' should be a struct with fields 'x' and 'names'
            % as output by modelspec::group_design
            write_covariate_names(design.names, outdir);
            write_matrix(design.x, fullfile(outdir, 'covar.dat'));
            % ^^ fnirs1 private functions
        end
    end
end
