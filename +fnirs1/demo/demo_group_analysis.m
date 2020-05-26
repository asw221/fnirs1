
%% Example group analysis 
% Group analyses require multiple participant data files and some guiding
% information about group level covariates. You can specify this group
% level information using data tables and an associated model formula
%

%% Locate and organize participant data files
% Retrieve data files from multiple directories and read patient
% demographic data
data_files = fnirs1.get_data_files;
data_files = [data_files; fnirs1.get_data_files];
data_files = reshape(data_files, numel(data_files), 1);

demo = readtable('~/Box/BayesianDataAnalysis/EN_CH_MA_PA_Tasks_NEW/Demographicinfo.xlsx');

% Subset the data and add a categorical column 
demo = demo(demo.ID > 3000, :);
demo = fnirs1.expand_table_conditions(demo, 2, 'Task');

% Check files are in the correct order:
[ok, mismatch] = fnirs1.check_filenames_id(data_files, demo.ID);


%% Fit the model
% Demographic information: names Cond and TempDeriv are special/reserved
% and do not need to be part of table. In a model formula, put Cond last of
% all categorical variables: e.g. Task + Cond + Task:Cond
%

fit = fnirs1.dlm(data_files, ...
    'GroupData', demo, ...
    'GroupFormula', 'ID ~ Task * Cond + LWIDraw + Age', ...
    'DownSampleRate', 10, ...
    'SpecificChannels', 9, ...
    'McmcControl', fnirs1.mcmc_control(30, false));

