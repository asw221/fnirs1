
data_files = fnirs1.get_data_files;
data_files = [data_files, fnirs1.get_data_files];
data_files = reshape(data_files, numel(data_files), 1);

demo = readtable('~/Box/BayesianDataAnalysis/EN_CH_MA_PA_Tasks_NEW/Demographicinfo.xlsx');

demo = demo(demo.ID > 3000, :);
demo = fnirs1.expand_table_conditions(demo, 2, 'Task');

% Check files are in the correct order:
[ok, mismatch] = fnirs1.check_filenames_id(data_files, demo.ID);

% demo = fnirs1.expand_table_conditions(demo, 3, 'Condition');


% Demographic information: names Cond and TempDeriv are special/reserved
% and do not need to be part of table. In model formula, put Cond last of
% all categorical variables: Task + Cond + Task:Cond
setupdir = fnirs1.specify_model(data_files, ...
    'GroupData', demo, ...
    'GroupFormula', 'ID ~ Task * Cond + LWIDraw + Age', ...
    'DownSampleRate', 10, ...
    'SpecificChannels', 9, ...
    'McmcControl', fnirs1.mcmc_debug);

setups = fnirs1.list_setup_files(setupdir);

% Group model - modular construction of setup file(s)
fit = fnirs1.dlm(setups);


% Or single subject - all at once
fit = fnirs1.dlm('', 'DownSampleRate', 10, 'SpecificChannels', 9, ...
    'McmcControl', fnirs1.mcmc_debug);

