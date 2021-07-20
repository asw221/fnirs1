
%% Create list of participant data files
% When multiple files appear in a row, fnirs1 will infer that these files
% correspond to subject replicates---the only thing you need to do to
% specify a repeated measures design like this is to organize your data
% files into this matrix-like format. Please note that fnirs1 does not
% support unbalanced designs, so each patient must have the same number of
% replicates
%
datadir = '~/Box/BayesianDataAnalysis/ConvertedData_20210609';

subjects = { ...
    'MA/CH1/1019_EN_MA_CH1.mat',    'PA/CH1/1019_EN_PA_CH1.mat'; ...
    'MA/CH1/1020_EN_MA_CH1.mat',    'PA/CH1/1020_EN_PA_CH1.mat'; ...
    'MA/CH1/1021_EN_MA_CH1.mat',    'PA/CH1/1021_EN_PA_CH1.mat'; ...
    'MA/CH1/1022_EN_MA_CH1.mat',    'PA/CH1/1022_EN_PA_CH1.mat'; ...
    'MA/CH1/1023_EN_MA_CH1.mat',    'PA/CH1/1023_EN_PA_CH1.mat'; ...
    'MA/CH1/1024_EN_MA_CH1.mat',    'PA/CH1/1024_EN_PA_CH1.mat'; ...
    'MA/CH1/1025_EN_MA_CH1.mat',    'PA/CH1/1025_EN_PA_CH1.mat'; ...
    'MA/CH1/1026_EN_MA_CH1.mat',    'PA/CH1/1026_EN_PA_CH1.mat'; ...
    'MA/CH1/1027_EN_MA_CH1.mat',    'PA/CH1/1027_EN_PA_CH1.mat'; ...
    'MA/CH1/1028_EN_MA_CH1.mat',    'PA/CH1/1028_EN_PA_CH1.mat'; ...
    'MA/CH1/1029_EN_MA_CH1.mat',    'PA/CH1/1029_EN_PA_CH1.mat'; ...
    'MA/CH1/1030_EN_MA_CH1.mat',    'PA/CH1/1030_EN_PA_CH1.mat'; ...
    'MA/CH1/1032_EN_MA_CH1.mat',    'PA/CH1/1032_EN_PA_CH1.mat'; ...
    'MA/CH1/1033_EN_MA_CH1.mat',    'PA/CH1/1033_EN_PA_CH1.mat'; ...
    'MA/CH1/1034_EN_MA_CH1.mat',    'PA/CH1/1034_EN_PA_CH1.mat' ...
    };

% Append full file path to each file in subjects
for i = 1:numel(subjects)
    subjects{i} = fullfile(datadir, subjects{i});
end


% Create demographic information table with task condition label. You
% likely already have demographic info in a table you can just read in
n = size(subjects, 1);
id = reshape([(1:n)', (1:n)']', [2*n, 1]);
pa = reshape([zeros([n 1]), ones([n, 1])]', [2*n, 1]);
demog = table(id, pa, 'VariableNames', {'ID', 'PAtask'});
% Note: dimension of demog is (2*n x 2): there are 2*n rows since each
% participant has 2 replicates. Variables like age, sex, etc. that are
% constant within patient should be repeated for subject replicate rows.
% See fnirs1.expand_table_conditions for a convenience method to repeat
% rows of a table if needed
%


%% Fit model using fnirs1's modular model setup
%

spec = fnirs1.specify_model(subjects, ...
    'DownSampleRate', 5, ...
    'McmcControl', fnirs1.mcmc_control(100), ...
    'GroupData', demog, ...
    'GroupFormula', 'ID ~ PAtask * Cond');
% Just as before, you can pass these arguments directly to fnirs1.dlm. In a
% small change from before, now you MUST supply 'GroupData' and
% 'GroupFormula' parameter pairs for group level analyses, or fnirs1 will
% complain at you


fit = fnirs1.dlm(spec);


%% Posterior credible intervals
%

% Add some contrasts:
fit = fit.add_contrast([-1, 0, 1]);

credint(fit, [0.8, 0.9, 0.95, 0.99])  % 80%, ..., 99% credible intervals

% To select specific channel (default is always the 1st channel in your
% subset):
% credint(fit, [0.8, 0.95], 10)  % CIs for 10th channel in subset
%
% Or for more decimal places:
credint(fit, [0.8, 0.95], 1, 2)  % CIs for 1st channel shown to 2 decimals


