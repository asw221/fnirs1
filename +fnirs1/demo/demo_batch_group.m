
%%             Example group analysis for batch execution
% Most of fnirs1 functionality is available directly from the command line.
% The only exceptions are that you won't have access to functions that
% explicitly rely on the Matlab GUI, like the file selector tool,
% fnirs1.get_data_files. This puts a slightly higher burden on the user to
% make sure data files names are correct, and any demographic variables are
% adequately wrangled and correctly ordered, etc.
%

% -> Edit the paths to the different data files to suit your own experiment

output_file_name = '/path/to/save/output.mat';

group_demographics_file = '/path/to/group/demographics/file.csv';
demographics = readtable(group_demographics_file);

% -> Reorganize demographics table if necessary

% Preprocessed fnirs data files with saved structures containing fields
% {hbo, hbr, hbt, ml, s, t} for outcomes, channel mappings, task stimulus
% design, and timepoints
fnirs_mat_files = { ...
    '/path/to/participant1.mat', ...
    '/path/to/participant2.mat', ...
    '/path/to/participant3.mat' ...
    };

% If the fnirs_mat_files' names contain the participants' study IDs and
% group_demographics_file contains a field with the same substring, you can
% uncomment the below to force a check of subject/data alignment
%
% if ~fnirs1.check_filenames_id(fnirs_mat_files, demographics.ID)
%     fnirs_mat_files = reshape(fnirs_mat_files, numel(fnirs_mat_files), 1);
%     [~, mm] = fnirs1.check_filenames_id(fnirs_mat_files, demographics.ID);
%     error('Mismatched files:\n%s', ...
%         sprintf('\t%s  -  %s\n', ...
%         [string(fnirs1.utils.basename(fnirs_mat_files(mm))), ...
%         string(demo.ID(mm))]') ...
%         );
% end


% -> Adjust model particulars

fit = fnirs1.dlm(fnirs_mat_files, ...
    'GroupData', fnirs1.utils.center(demographics), ...
    'GroupFormula', 'ID ~ LWIDraw + Age + Cond', ...
    'DownSampleRate', 5, ...
    'SpecificChannels', 9:10);

disp(geweke(fit));  % Print rough diagnostic information

% -> Add contrasts, etc if desired


save(output_file_name, 'fit', 'fnirs_mat_files', 'group_demographics_file');


