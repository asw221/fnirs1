
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
        Intervals;     % Posterior credible intervals for each parameter
        Level;         % Credible interval probability width
        StdErrors;     % Standard error of each Estimate
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
    methods
        function obj = dlm_summary()
            % Default constructor for FNIRS1.DLM_SUMMARY objects
            obj.Nickname = '';
            obj.Descriptions = {};
            obj.Estimates = [];
            obj.Intervals = [];
            obj.Level = 0.95;
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
        end
        function obj = read_from_file(obj, file)
            % Read data from a Parameter_Estimates.log file and extract
            % parameter estimates and inferential summaries
            
            % Read file into cell array and clean up lines
            lines = fnirs1.read_whole_file(file);
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
                        sprintf('Cond %i', condNo + 1));
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
            end
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
                        descrip = sprintf(obj(j).FormatDescrip, obj(j).Descriptions{i});
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
    end
    methods (Static)
        function displayEmptyObject()
            fprintf('\tempty fnirs1.dlm_summary object\n');
        end
    end
end
